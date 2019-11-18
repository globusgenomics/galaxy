from __future__ import absolute_import

import base64
import json
import logging

from markupsafe import escape
from six.moves.html_parser import HTMLParser
from six.moves.http_client import HTTPConnection
from sqlalchemy import and_
from sqlalchemy.orm import eagerload, joinedload, lazyload, undefer
from sqlalchemy.sql import expression

from galaxy import (
    model,
    util,
    web
)
from galaxy.managers import workflows
from galaxy.model.item_attrs import UsesItemRatings
from galaxy.model.mapping import desc
from galaxy.security.validate_user_input import validate_publicname
from galaxy.tools.parameters.basic import workflow_building_modes, RuntimeValue, DataToolParameter
from galaxy.util import (
    FILENAME_VALID_CHARS,
    unicodify
)
from galaxy.util.sanitize_html import sanitize_html
from galaxy.web import error, url_for
from galaxy.web.base.controller import (
    BaseUIController,
    SharableMixin,
    UsesStoredWorkflowMixin
)
from galaxy.web.framework.helpers import (
    grids,
    time_ago,
)
from galaxy.workflow.extract import (
    extract_workflow,
    summarize
)
from galaxy.workflow.modules import (
    module_factory,
    WorkflowModuleInjector
)
from galaxy.workflow.render import (
    STANDALONE_SVG_TEMPLATE,
    WorkflowCanvas
)

log = logging.getLogger(__name__)


class StoredWorkflowListGrid(grids.Grid):

    class StepsColumn(grids.GridColumn):
        def get_value(self, trans, grid, workflow):
            return len(workflow.latest_workflow.steps)

    # Grid definition
    use_panels = True
    title = "Saved Workflows"
    model_class = model.StoredWorkflow
    default_filter = {"name": "All", "tags": "All"}
    default_sort_key = "-update_time"
    columns = [
        grids.TextColumn("Name", key="name", attach_popup=True, filterable="advanced"),
        grids.IndividualTagsColumn("Tags",
                                   "tags",
                                   model_tag_association_class=model.StoredWorkflowTagAssociation,
                                   filterable="advanced",
                                   grid_name="StoredWorkflowListGrid"),
        StepsColumn("Steps"),
        grids.GridColumn("Created", key="create_time", format=time_ago),
        grids.GridColumn("Last Updated", key="update_time", format=time_ago),
    ]
    columns.append(
        grids.MulticolFilterColumn(
            "Search",
            cols_to_filter=[columns[0], columns[1]],
            key="free-text-search", visible=False, filterable="standard"
        )
    )
    operations = [
        grids.GridOperation("Edit", allow_multiple=False, condition=(lambda item: not item.deleted), async_compatible=False),
        grids.GridOperation("Run", condition=(lambda item: not item.deleted), async_compatible=False),
        grids.GridOperation("Copy", condition=(lambda item: not item.deleted), async_compatible=False),
        grids.GridOperation("Rename", condition=(lambda item: not item.deleted), async_compatible=False),
        grids.GridOperation("Sharing", condition=(lambda item: not item.deleted), async_compatible=False),
        grids.GridOperation("Delete", condition=(lambda item: item.deleted), async_compatible=True),
    ]

    def apply_query_filter(self, trans, query, **kwargs):
        return query.filter_by(user=trans.user, deleted=False)


class StoredWorkflowAllPublishedGrid(grids.Grid):
    title = "Published Workflows"
    model_class = model.StoredWorkflow
    default_sort_key = "update_time"
    default_filter = dict(public_url="All", username="All", tags="All")
    columns = [
        grids.PublicURLColumn("Name", key="name", filterable="advanced", attach_popup=True),
        grids.OwnerAnnotationColumn("Annotation",
                                    key="annotation",
                                    model_annotation_association_class=model.StoredWorkflowAnnotationAssociation,
                                    filterable="advanced"),
        grids.OwnerColumn("Owner", key="username", model_class=model.User, filterable="advanced"),
        grids.CommunityRatingColumn("Community Rating", key="rating"),
        grids.CommunityTagsColumn("Community Tags", key="tags",
                                  model_tag_association_class=model.StoredWorkflowTagAssociation,
                                  filterable="advanced", grid_name="PublicWorkflowListGrid"),
        grids.ReverseSortColumn("Last Updated", key="update_time", format=time_ago)
    ]
    columns.append(
        grids.MulticolFilterColumn(
            "Search name, annotation, owner, and tags",
            cols_to_filter=[columns[0], columns[1], columns[2], columns[4]],
            key="free-text-search", visible=False, filterable="standard"
        )
    )
    operations = [
        grids.GridOperation(
            'Run',
            condition=(lambda item: not item.deleted),
            allow_multiple=False,
            url_args=dict(controller='workflows', action="run")
        ),
        grids.GridOperation(
            "Import",
            condition=(lambda item: not item.deleted),
            allow_multiple=False,
            url_args=dict(action="imp")
        ),
        grids.GridOperation(
            "Save as File",
            condition=(lambda item: not item.deleted),
            allow_multiple=False,
            url_args=dict(action="export_to_file")
        ),
    ]
    num_rows_per_page = 50
    use_paging = True

    def build_initial_query(self, trans, **kwargs):
        # See optimization description comments and TODO for tags in matching public histories query.
        # In addition to that - be sure to lazyload the latest_workflow - it isn't needed and it causes all
        # of its steps to be eagerly loaded.
        return trans.sa_session.query(self.model_class).join("user").options(lazyload("latest_workflow"), eagerload("user").load_only("username"), eagerload("annotations"), undefer("average_rating"))

    def apply_query_filter(self, trans, query, **kwargs):
        # A public workflow is published, has a slug, and is not deleted.
        return query.filter(
            self.model_class.published == expression.true()).filter(
            self.model_class.slug.isnot(None)).filter(
            self.model_class.deleted == expression.false())


# Simple HTML parser to get all content in a single tag.
class SingleTagContentsParser(HTMLParser):

    def __init__(self, target_tag):
        # Cannot use super() because HTMLParser is an old-style class in Python2
        HTMLParser.__init__(self)
        self.target_tag = target_tag
        self.cur_tag = None
        self.tag_content = ""

    def handle_starttag(self, tag, attrs):
        """ Called for each start tag. """
        self.cur_tag = tag

    def handle_data(self, text):
        """ Called for each block of plain text. """
        if self.cur_tag == self.target_tag:
            self.tag_content += text


class WorkflowController(BaseUIController, SharableMixin, UsesStoredWorkflowMixin, UsesItemRatings):
    stored_list_grid = StoredWorkflowListGrid()
    published_list_grid = StoredWorkflowAllPublishedGrid()

    __myexp_url = "www.myexperiment.org:80"

    @web.expose
    @web.require_login("use Galaxy workflows")
    def list_grid(self, trans, **kwargs):
        """ List user's stored workflows. """
        # status = message = None
        if 'operation' in kwargs:
            operation = kwargs['operation'].lower()
            if operation == "rename":
                return self.rename(trans, **kwargs)
            history_ids = util.listify(kwargs.get('id', []))
            if operation == "sharing":
                return self.sharing(trans, id=history_ids)
        return self.stored_list_grid(trans, **kwargs)

    @web.expose
    @web.require_login("use Galaxy workflows", use_panels=True)
    def list(self, trans):
        """
        Render workflow main page (management of existing workflows)
        """
        # Take care of proxy prefix in url as well
        redirect_url = url_for('/') + 'workflow'
        return trans.response.send_redirect(redirect_url)

    @web.expose
    @web.json
    def list_published(self, trans, **kwargs):
        return self.published_list_grid(trans, **kwargs)

    @web.expose
    def display_by_username_and_slug(self, trans, username, slug, format='html'):
        """
        Display workflow based on a username and slug. Format can be html, json, or json-download.
        """

        # Get workflow by username and slug. Security is handled by the display methods below.
        session = trans.sa_session
        user = session.query(model.User).filter_by(username=username).first()
        if not user:
            raise web.httpexceptions.HTTPNotFound()
        stored_workflow = trans.sa_session.query(model.StoredWorkflow).filter_by(user=user, slug=slug, deleted=False).first()
        if not stored_workflow:
            raise web.httpexceptions.HTTPNotFound()
        encoded_id = trans.security.encode_id(stored_workflow.id)

        # Display workflow in requested format.
        if format == 'html':
            return self._display(trans, stored_workflow)
        elif format == 'json':
            return self.for_direct_import(trans, encoded_id)
        elif format == 'json-download':
            return self.export_to_file(trans, encoded_id)
        elif format == 'export-parameters':
            return self.export_parameters( trans, encoded_id )

    @web.expose
    def display_by_id(self, trans, id):
        """ Display workflow based on id. """
        # Get workflow.
        stored_workflow = self.get_stored_workflow(trans, id)
        return self._display(trans, stored_workflow)

    def _display(self, trans, stored_workflow):
        """ Diplay workflow as HTML page. """

        if stored_workflow is None:
            raise web.httpexceptions.HTTPNotFound()
        # Security check raises error if user cannot access workflow.
        self.security_check(trans, stored_workflow, False, True)
        # Get data for workflow's steps.
        self.get_stored_workflow_steps(trans, stored_workflow)
        # Get annotations.
        stored_workflow.annotation = self.get_item_annotation_str(trans.sa_session, stored_workflow.user, stored_workflow)
        for step in stored_workflow.latest_workflow.steps:
            step.annotation = self.get_item_annotation_str(trans.sa_session, stored_workflow.user, step)
        # Get rating data.
        user_item_rating = 0
        if trans.get_user():
            user_item_rating = self.get_user_item_rating(trans.sa_session, trans.get_user(), stored_workflow)
            if user_item_rating:
                user_item_rating = user_item_rating.rating
            else:
                user_item_rating = 0
        ave_item_rating, num_ratings = self.get_ave_item_rating_data(trans.sa_session, stored_workflow)
        return trans.fill_template_mako("workflow/display.mako", item=stored_workflow, item_data=stored_workflow.latest_workflow.steps,
                                        user_item_rating=user_item_rating, ave_item_rating=ave_item_rating, num_ratings=num_ratings)

    @web.expose
    def get_item_content_async(self, trans, id):
        """ Returns item content in HTML format. """

        stored = self.get_stored_workflow(trans, id, False, True)
        if stored is None:
            raise web.httpexceptions.HTTPNotFound()

        # Get data for workflow's steps.
        self.get_stored_workflow_steps(trans, stored)
        # Get annotations.
        stored.annotation = self.get_item_annotation_str(trans.sa_session, stored.user, stored)
        for step in stored.latest_workflow.steps:
            step.annotation = self.get_item_annotation_str(trans.sa_session, stored.user, step)
        return trans.stream_template_mako("/workflow/item_content.mako", item=stored, item_data=stored.latest_workflow.steps)

    @web.expose
    @web.require_login("use Galaxy workflows")
    def share(self, trans, id, email="", use_panels=False):
        msg = mtype = None
        # Load workflow from database
        stored = self.get_stored_workflow(trans, id)
        if email:
            other = trans.sa_session.query(model.User) \
                                    .filter(and_(model.User.table.c.email == email,
                                                 model.User.table.c.deleted == expression.false())) \
                                    .first()
            if not other:
                mtype = "error"
                msg = ("User '%s' does not exist" % escape(email))
            elif other == trans.get_user():
                mtype = "error"
                msg = ("You cannot share a workflow with yourself")
            elif trans.sa_session.query(model.StoredWorkflowUserShareAssociation) \
                    .filter_by(user=other, stored_workflow=stored).count() > 0:
                mtype = "error"
                msg = ("Workflow already shared with '%s'" % escape(email))
            else:
                share = model.StoredWorkflowUserShareAssociation()
                share.stored_workflow = stored
                share.user = other
                session = trans.sa_session
                session.add(share)
                session.flush()
                trans.set_message("Workflow '%s' shared with user '%s'" % (escape(stored.name), escape(other.email)))
                return trans.response.send_redirect(url_for(controller='workflow', action='sharing', id=id))
        return trans.fill_template("/ind_share_base.mako",
                                   message=msg,
                                   messagetype=mtype,
                                   item=stored,
                                   email=email,
                                   use_panels=use_panels)

    @web.expose
    @web.require_login("Share or export Galaxy workflows")
    def sharing(self, trans, id, **kwargs):
        """ Handle workflow sharing. """
        session = trans.sa_session
        if 'unshare_me' in kwargs:
            # Remove self from shared associations with workflow.
            stored = self.get_stored_workflow(trans, id, False, True)
            association = session.query(model.StoredWorkflowUserShareAssociation) \
                                 .filter_by(user=trans.user, stored_workflow=stored).one()
            session.delete(association)
            session.flush()
            return self.list(trans)
        else:
            # Get session and workflow.
            stored = self.get_stored_workflow(trans, id)
            session.add(stored)

            # Do operation on workflow.
            if 'make_accessible_via_link' in kwargs:
                self._make_item_accessible(trans.sa_session, stored)
            elif 'make_accessible_and_publish' in kwargs:
                self._make_item_accessible(trans.sa_session, stored)
                stored.published = True
            elif 'publish' in kwargs:
                stored.published = True
            elif 'disable_link_access' in kwargs:
                stored.importable = False
            elif 'unpublish' in kwargs:
                stored.published = False
            elif 'disable_link_access_and_unpublish' in kwargs:
                stored.importable = stored.published = False
            elif 'unshare_user' in kwargs:
                user = session.query(model.User).get(trans.security.decode_id(kwargs['unshare_user']))
                if not user:
                    error("User not found for provided id")
                association = session.query(model.StoredWorkflowUserShareAssociation) \
                                     .filter_by(user=user, stored_workflow=stored).one()
                session.delete(association)

            # Legacy issue: workflows made accessible before recent updates may not have a slug. Create slug for any workflows that need them.
            if stored.importable and not stored.slug:
                self._make_item_accessible(trans.sa_session, stored)

            session.flush()
            return trans.fill_template("/workflow/sharing.mako", use_panels=True, item=stored)

    @web.expose
    @web.require_login("share Galaxy items")
    def set_public_username(self, trans, id, username, **kwargs):
        """ Set user's public username and delegate to sharing() """
        user = trans.get_user()
        # message from validate_publicname does not contain input, no need
        # to escape.
        message = validate_publicname(trans, username, user)
        if message:
            return trans.fill_template("/workflow/sharing.mako", item=self.get_item(trans, id), message=message, status="error")
        user.username = username
        trans.sa_session.flush()
        return self.sharing(trans, id, **kwargs)

    @web.expose
    @web.require_login("to import a workflow", use_panels=True)
    def imp(self, trans, id, **kwargs):
        """Imports a workflow shared by other users."""
        # Set referer message.
        referer = trans.request.referer
        if referer:
            referer_message = "<a href='%s'>return to the previous page</a>" % escape(referer)
        else:
            referer_message = "<a href='%s'>go to Galaxy's start page</a>" % url_for('/')

        # Do import.
        stored = self.get_stored_workflow(trans, id, check_ownership=False)
        if stored.importable is False:
            return trans.show_error_message("The owner of this workflow has disabled imports via this link.<br>You can %s" % referer_message, use_panels=True)
        elif stored.deleted:
            return trans.show_error_message("You can't import this workflow because it has been deleted.<br>You can %s" % referer_message, use_panels=True)
        self._import_shared_workflow(trans, stored)

        # Redirect to load galaxy frames.
        return trans.show_ok_message(
            message="""Workflow "%s" has been imported. <br>You can <a href="%s">start using this workflow</a> or %s."""
            % (stored.name, web.url_for('/workflows/list'), referer_message))

    @web.expose
    @web.require_login("use Galaxy workflows")
    def rename_async(self, trans, id, new_name=None, **kwargs):
        stored = self.get_stored_workflow(trans, id)
        if new_name:
            san_new_name = sanitize_html(new_name)
            stored.name = san_new_name
            stored.latest_workflow.name = san_new_name
            trans.sa_session.flush()
            return stored.name

    @web.expose
    @web.require_login("use Galaxy workflows")
    def annotate_async(self, trans, id, new_annotation=None, **kwargs):
        stored = self.get_stored_workflow(trans, id)
        if new_annotation:
            # Sanitize annotation before adding it.
            new_annotation = sanitize_html(new_annotation)
            self.add_item_annotation(trans.sa_session, trans.get_user(), stored, new_annotation)
            trans.sa_session.flush()
            return new_annotation

    @web.expose
    @web.require_login("rate items")
    @web.json
    def rate_async(self, trans, id, rating):
        """ Rate a workflow asynchronously and return updated community data. """

        stored = self.get_stored_workflow(trans, id, check_ownership=False, check_accessible=True)
        if not stored:
            return trans.show_error_message("The specified workflow does not exist.")

        # Rate workflow.
        self.rate_item(trans.sa_session, trans.get_user(), stored, rating)

        return self.get_ave_item_rating_data(trans.sa_session, stored)

    @web.expose
    @web.require_login("use Galaxy workflows")
    def set_accessible_async(self, trans, id=None, accessible=False):
        """ Set workflow's importable attribute and slug. """
        stored = self.get_stored_workflow(trans, id)

        # Only set if importable value would change; this prevents a change in the update_time unless attribute really changed.
        importable = accessible in ['True', 'true', 't', 'T']
        if stored and stored.importable != importable:
            if importable:
                self._make_item_accessible(trans.sa_session, stored)
            else:
                stored.importable = importable
            trans.sa_session.flush()
        return

    @web.expose
    def get_embed_html_async(self, trans, id):
        """ Returns HTML for embedding a workflow in a page. """

        # TODO: user should be able to embed any item he has access to. see display_by_username_and_slug for security code.
        stored = self.get_stored_workflow(trans, id)
        if stored:
            return "Embedded Workflow '%s'" % stored.name

    @web.expose
    @web.json
    @web.require_login("use Galaxy workflows")
    def get_name_and_link_async(self, trans, id=None):
        """ Returns workflow's name and link. """
        stored = self.get_stored_workflow(trans, id)

        return_dict = {"name": stored.name,
                       "link": url_for(controller='workflow',
                                       action="display_by_username_and_slug",
                                       username=stored.user.username,
                                       slug=stored.slug)}
        return return_dict

    @web.expose
    @web.require_login("use Galaxy workflows")
    def gen_image(self, trans, id):
        stored = self.get_stored_workflow(trans, id, check_ownership=True)
        try:
            svg = self._workflow_to_svg_canvas(trans, stored)
        except Exception:
            status = 'error'
            message = 'Galaxy is unable to create the SVG image. Please check your workflow, there might be missing tools.'
            return trans.fill_template("/workflow/sharing.mako", use_panels=True, item=stored, status=status, message=message)
        trans.response.set_content_type("image/svg+xml")
        s = STANDALONE_SVG_TEMPLATE % svg.tostring()
        return s.encode('utf-8')

    @web.expose
    @web.require_login("use Galaxy workflows")
    def copy(self, trans, id, save_as_name=None):
        # Get workflow to copy.
        stored = self.get_stored_workflow(trans, id, check_ownership=False)
        user = trans.get_user()
        if stored.user == user:
            owner = True
        else:
            if trans.sa_session.query(model.StoredWorkflowUserShareAssociation) \
                    .filter_by(user=user, stored_workflow=stored).count() == 0:
                error("Workflow is not owned by or shared with current user")
            owner = False

        # Copy.
        new_stored = model.StoredWorkflow()
        if (save_as_name):
            new_stored.name = '%s' % save_as_name
        else:
            new_stored.name = "Copy of %s" % stored.name
        new_stored.latest_workflow = stored.latest_workflow
        # Copy annotation.
        annotation_obj = self.get_item_annotation_obj(trans.sa_session, stored.user, stored)
        if annotation_obj:
            self.add_item_annotation(trans.sa_session, trans.get_user(), new_stored, annotation_obj.annotation)
        new_stored.copy_tags_from(trans.user, stored)
        if not owner:
            new_stored.name += " shared by %s" % stored.user.email
        new_stored.user = user
        # Persist
        session = trans.sa_session
        session.add(new_stored)
        session.flush()
        # Display the management page
        message = 'Created new workflow with name: %s' % escape(new_stored.name)
        trans.set_message(message)
        return_url = url_for('/') + 'workflow?status=done&message=%s' % escape(message)
        trans.response.send_redirect(return_url)

    @web.legacy_expose_api
    def create(self, trans, payload=None, **kwd):
        if trans.request.method == 'GET':
            return {
                'title'  : 'Create Workflow',
                'inputs' : [{
                    'name'  : 'workflow_name',
                    'label' : 'Name',
                    'value' : 'Unnamed workflow'
                }, {
                    'name'  : 'workflow_annotation',
                    'label' : 'Annotation',
                    'help'  : 'A description of the workflow; annotation is shown alongside shared or published workflows.'
                }]}
        else:
            user = trans.get_user()
            workflow_name = payload.get('workflow_name')
            workflow_annotation = payload.get('workflow_annotation')
            if not workflow_name:
                return self.message_exception(trans, 'Please provide a workflow name.')
            # Create the new stored workflow
            stored_workflow = model.StoredWorkflow()
            stored_workflow.name = workflow_name
            stored_workflow.user = user
            self.create_item_slug(trans.sa_session, stored_workflow)
            # And the first (empty) workflow revision
            workflow = model.Workflow()
            workflow.name = workflow_name
            workflow.stored_workflow = stored_workflow
            stored_workflow.latest_workflow = workflow
            # Add annotation.
            workflow_annotation = sanitize_html(workflow_annotation)
            self.add_item_annotation(trans.sa_session, trans.get_user(), stored_workflow, workflow_annotation)
            # Persist
            session = trans.sa_session
            session.add(stored_workflow)
            session.flush()
            return {'id': trans.security.encode_id(stored_workflow.id), 'message': 'Workflow %s has been created.' % workflow_name}

    @web.json
    def save_workflow_as(self, trans, workflow_name, workflow_data, workflow_annotation=""):
        """
            Creates a new workflow based on Save As command. It is a new workflow, but
            is created with workflow_data already present.
        """
        user = trans.get_user()
        if workflow_name is not None:
            workflow_contents_manager = workflows.WorkflowContentsManager(trans.app)
            stored_workflow = model.StoredWorkflow()
            stored_workflow.name = workflow_name
            stored_workflow.user = user
            self.create_item_slug(trans.sa_session, stored_workflow)
            workflow = model.Workflow()
            workflow.name = workflow_name
            workflow.stored_workflow = stored_workflow
            stored_workflow.latest_workflow = workflow
            # Add annotation.
            workflow_annotation = sanitize_html(workflow_annotation)
            self.add_item_annotation(trans.sa_session, trans.get_user(), stored_workflow, workflow_annotation)

            # Persist
            session = trans.sa_session
            session.add(stored_workflow)
            session.flush()

            try:
                workflow, errors = workflow_contents_manager.update_workflow_from_raw_description(
                    trans,
                    stored_workflow,
                    workflow_data,
                )
            except workflows.MissingToolsException as e:
                return dict(
                    name=e.workflow.name,
                    message=("This workflow includes missing or invalid tools. "
                             "It cannot be saved until the following steps are removed or the missing tools are enabled."),
                    errors=e.errors,
                )
            return (trans.security.encode_id(stored_workflow.id))
        else:
            # This is an error state, 'save as' must have a workflow_name
            log.exception("Error in Save As workflow: no name.")

    @web.expose
    def delete(self, trans, id=None):
        """
        Mark a workflow as deleted
        """
        # Load workflow from database
        stored = self.get_stored_workflow(trans, id)
        # Mark as deleted and save
        stored.deleted = True
        trans.user.stored_workflow_menu_entries = [entry for entry in trans.user.stored_workflow_menu_entries if entry.stored_workflow != stored]
        trans.sa_session.add(stored)
        trans.sa_session.flush()
        # Display the management page
        message = "Workflow deleted: %s" % escape(stored.name)
        trans.set_message(message)
        return trans.response.send_redirect(url_for('/') + 'workflow?status=done&message=%s' % escape(message))

    @web.expose
    @web.require_login("edit workflows")
    def editor(self, trans, id=None, version=None):
        """
        Render the main workflow editor interface. The canvas is embedded as
        an iframe (necessary for scrolling to work properly), which is
        rendered by `editor_canvas`.
        """
        if not id:
            error("Invalid workflow id")
        stored = self.get_stored_workflow(trans, id)
        # The following query loads all user-owned workflows,
        # So that they can be copied or inserted in the workflow editor.
        workflows = trans.sa_session.query(model.StoredWorkflow) \
            .filter_by(user=trans.user, deleted=False) \
            .order_by(desc(model.StoredWorkflow.table.c.update_time)) \
            .options(joinedload('latest_workflow').joinedload('steps')) \
            .all()
        if version is None:
            version = len(stored.workflows) - 1
        else:
            version = int(version)
        return trans.fill_template("workflow/editor.mako",
                                   workflows=workflows,
                                   stored=stored,
                                   version=version,
                                   annotation=self.get_item_annotation_str(trans.sa_session, trans.user, stored))

    @web.json
    def load_workflow(self, trans, id, version=None):
        """
        Get the latest Workflow for the StoredWorkflow identified by `id` and
        encode it as a json string that can be read by the workflow editor
        web interface.
        """
        trans.workflow_building_mode = workflow_building_modes.ENABLED
        stored = self.get_stored_workflow(trans, id, check_ownership=True, check_accessible=False)
        workflow_contents_manager = workflows.WorkflowContentsManager(trans.app)
        return workflow_contents_manager.workflow_to_dict(trans, stored, style="editor", version=version)

    @web.expose
    @web.require_login("use workflows")
    def export_to_myexp(self, trans, id, myexp_username, myexp_password):
        """
        Exports a workflow to myExperiment website.
        """
        trans.workflow_building_mode = workflow_building_modes.ENABLED
        stored = self.get_stored_workflow(trans, id, check_ownership=False, check_accessible=True)

        # Convert workflow to dict.
        workflow_dict = self._workflow_to_dict(trans, stored)

        #
        # Create and submit workflow myExperiment request.
        #

        # Create workflow content JSON.
        workflow_content = json.dumps(workflow_dict, indent=4, sort_keys=True)

        # Create myExperiment request.
        request_raw = trans.fill_template(
            "workflow/myexp_export.mako",
            workflow_name=workflow_dict['name'],
            workflow_description=workflow_dict['annotation'],
            workflow_content=workflow_content,
            workflow_svg=self._workflow_to_svg_canvas(trans, stored).tostring()
        )
        # strip() b/c myExperiment XML parser doesn't allow white space before XML; utf-8 handles unicode characters.
        request = unicodify(request_raw.strip(), 'utf-8')

        # Do request and get result.
        auth_header = base64.b64encode('%s:%s' % (myexp_username, myexp_password))
        headers = {"Content-type": "text/xml", "Accept": "text/xml", "Authorization": "Basic %s" % auth_header}
        myexp_url = trans.app.config.get("myexperiment_url", self.__myexp_url)
        conn = HTTPConnection(myexp_url)
        # NOTE: blocks web thread.
        conn.request("POST", "/workflow.xml", request, headers)
        response = conn.getresponse()
        response_data = response.read()
        conn.close()

        # Do simple parse of response to see if export successful and provide user feedback.
        parser = SingleTagContentsParser('id')
        parser.feed(response_data)
        myexp_workflow_id = parser.tag_content
        workflow_list_str = " <br>Return to <a href='%s'>workflow list." % url_for(controller='workflows', action='list')
        if myexp_workflow_id:
            return trans.show_message(
                """Workflow '%s' successfully exported to myExperiment. <br/>
                <a href="http://%s/workflows/%s">Click here to view the workflow on myExperiment</a> %s
                """ % (stored.name, myexp_url, myexp_workflow_id, workflow_list_str),
                use_panels=True)
        else:
            return trans.show_error_message(
                "Workflow '%s' could not be exported to myExperiment. Error: %s %s" %
                (stored.name, response_data, workflow_list_str), use_panels=True)

    @web.json_pretty
    def for_direct_import(self, trans, id):
        """
        Get the latest Workflow for the StoredWorkflow identified by `id` and
        encode it as a json string that can be imported back into Galaxy

        This has slightly different information than the above. In particular,
        it does not attempt to decode forms and build UIs, it just stores
        the raw state.
        """
        stored = self.get_stored_workflow(trans, id, check_ownership=False, check_accessible=True)
        return self._workflow_to_dict(trans, stored)

    @web.json_pretty
    def export_to_file(self, trans, id):
        """
        Get the latest Workflow for the StoredWorkflow identified by `id` and
        encode it as a json string that can be imported back into Galaxy

        This has slightly different information than the above. In particular,
        it does not attempt to decode forms and build UIs, it just stores
        the raw state.
        """

        # Get workflow.
        stored = self.get_stored_workflow(trans, id, check_ownership=False, check_accessible=True)

        # Stream workflow to file.
        stored_dict = self._workflow_to_dict(trans, stored)
        if not stored_dict:
            # This workflow has a tool that's missing from the distribution
            trans.response.status = 400
            return "Workflow cannot be exported due to missing tools."
        sname = stored.name
        sname = ''.join(c in FILENAME_VALID_CHARS and c or '_' for c in sname)[0:150]
        trans.response.headers["Content-Disposition"] = 'attachment; filename="Galaxy-Workflow-%s.ga"' % (sname)
        trans.response.set_content_type('application/galaxy-archive')
        return stored_dict

    @web.expose
    def build_from_current_history(self, trans, job_ids=None, dataset_ids=None, dataset_collection_ids=None, workflow_name=None, dataset_names=None, dataset_collection_names=None):
        user = trans.get_user()
        history = trans.get_history()
        if not user:
            return trans.show_error_message("Must be logged in to create workflows")
        if (job_ids is None and dataset_ids is None) or workflow_name is None:
            jobs, warnings = summarize(trans)
            # Render
            return trans.fill_template(
                "workflow/build_from_current_history.mako",
                jobs=jobs,
                warnings=warnings,
                history=history
            )
        else:
            # If there is just one dataset name selected or one dataset collection, these
            # come through as string types instead of lists. xref #3247.
            dataset_names = util.listify(dataset_names)
            dataset_collection_names = util.listify(dataset_collection_names)
            stored_workflow = extract_workflow(
                trans,
                user=user,
                job_ids=job_ids,
                dataset_ids=dataset_ids,
                dataset_collection_ids=dataset_collection_ids,
                workflow_name=workflow_name,
                dataset_names=dataset_names,
                dataset_collection_names=dataset_collection_names
            )
            # Index page with message
            workflow_id = trans.security.encode_id(stored_workflow.id)
            return trans.show_message('Workflow "%s" created from current history. '
                                      'You can <a href="%s" target="_parent">edit</a> or <a href="%s" target="_parent">run</a> the workflow.'
                                      % (escape(workflow_name), url_for(controller='workflow', action='editor', id=workflow_id),
                                         url_for(controller='workflows', action='run', id=workflow_id)))

    def get_item(self, trans, id):
        return self.get_stored_workflow(trans, id)

    @web.expose
    def tag_outputs(self, trans, id, **kwargs):
        stored = self.get_stored_workflow(trans, id, check_ownership=False)
        user = trans.get_user()
        if stored.user != user:
            if trans.sa_session.query(model.StoredWorkflowUserShareAssociation) \
                    .filter_by(user=user, stored_workflow=stored).count() == 0:
                error("Workflow is not owned by or shared with current user")
        # Get the latest revision
        workflow = stored.latest_workflow
        # It is possible for a workflow to have 0 steps
        if len(workflow.steps) == 0:
            error("Workflow cannot be tagged for outputs because it does not have any steps")
        if workflow.has_cycles:
            error("Workflow cannot be tagged for outputs because it contains cycles")
        if workflow.has_errors:
            error("Workflow cannot be tagged for outputs because of validation errors in some steps")
        # Build the state for each step
        errors = {}
        has_upgrade_messages = False
        # has_errors is never used
        # has_errors = False
        if kwargs:
            # If kwargs were provided, the states for each step should have
            # been POSTed
            for step in workflow.steps:
                if step.type == 'tool':
                    # Extract just the output flags for this step.
                    p = "%s|otag|" % step.id
                    l = len(p)
                    outputs = [k[l:] for (k, v) in kwargs.items() if k.startswith(p)]
                    if step.workflow_outputs:
                        for existing_output in step.workflow_outputs:
                            if existing_output.output_name not in outputs:
                                trans.sa_session.delete(existing_output)
                            else:
                                outputs.remove(existing_output.output_name)
                    for outputname in outputs:
                        m = model.WorkflowOutput(workflow_step_id=int(step.id), output_name=outputname)
                        trans.sa_session.add(m)
        # Prepare each step
        trans.sa_session.flush()
        module_injector = WorkflowModuleInjector(trans)
        for step in workflow.steps:
            step.upgrade_messages = {}
            # Contruct modules
            module_injector.inject(step)
            if step.upgrade_messages:
                has_upgrade_messages = True
            if step.type == 'tool' or step.type is None:
                # Error dict
                if step.tool_errors:
                    errors[step.id] = step.tool_errors
        # Render the form
        return trans.fill_template(
            "workflow/tag_outputs.mako",
            steps=workflow.steps,
            workflow=stored,
            has_upgrade_messages=has_upgrade_messages,
            errors=errors,
            incoming=kwargs
        )

    def _workflow_to_svg_canvas(self, trans, stored):
        workflow = stored.latest_workflow
        workflow_canvas = WorkflowCanvas()
        for step in workflow.steps:
            # Load from database representation
            module = module_factory.from_workflow_step(trans, step)
            module_name = module.get_name()
            module_data_inputs = module.get_data_inputs()
            module_data_outputs = module.get_data_outputs()
            workflow_canvas.populate_data_for_step(
                step,
                module_name,
                module_data_inputs,
                module_data_outputs,
            )
        workflow_canvas.add_steps()
        return workflow_canvas.finish()

    @web.expose
    @web.require_login( "use workflows" )
    def parameters( self, trans, id=None, **kwargs ):
        """ Handle workflow sharing. """
        stored = self.get_stored_workflow( trans, id, check_ownership=False, check_accessible=True )

        return trans.fill_template( "/workflow/export_parameters.mako", use_panels=True, item=stored )

    #@web.expose
    def export_parameters( self, trans, id, **kwargs  ):
        stored = self.get_stored_workflow( trans, id, check_ownership=False )
        user = trans.get_user()
        if stored.user != user:
            if trans.sa_session.query( model.StoredWorkflowUserShareAssociation ) \
                    .filter_by( user=user, stored_workflow=stored ).count() == 0:
                error( "Workflow is not owned by or shared with current user" )
        # Get the latest revision
        workflow = stored.latest_workflow
#        print workflow

#        for property, value in vars(workflow).iteritems():
#            print property, ": ", value

        # It is possible for a workflow to have 0 steps
        if len( workflow.steps ) == 0:
            error( "Workflow cannot be run because it does not have any steps" )
        #workflow = Workflow.from_simple( simplejson.loads( stored.encoded_value ), trans.app )
        if workflow.has_cycles:
            error( "Workflow cannot be run because it contains cycles" )
        if workflow.has_errors:
            error( "Workflow cannot be run because of validation errors in some steps" )
        # Build the state for each step
        errors = {}
        has_upgrade_messages = False
        has_errors = False
        saved_history = None
        data_export = ""
        data_export += "#Data Export for Workflow Batch Submission Through the APII\n\n\n"
        data_export += "### INSTRUCTIONS\n"
        data_export += "#######################################\n"
        data_export += "#The following data can be used to input the parameters you have previously determined to be\n"
        data_export += "#set at runtime. Please specify the library or history where the input data can be found.\n"
        data_export += "#Once you have filled out the table you can run the API script to submit the jobs through Galaxy\n"
        data_export += "#via the API.\n\n"
        data_export += "#NOTE: If you make any changes to the workflow or edit the name of the workflow, you will need\n"
        data_export += "#to recreate the table before submitting the job via the API since some metadata parameters will\n"
        data_export += "#be modified.\n\n"            
        data_export += "#NOTE: It is up to the user to make sure the input files are in the correct format for each\n"
        data_export += "#parameter being filled out.\n\n"
        data_export += "#NOTE: You will need to specify three items for input files to an application.\n"
        data_export += "#The format for an input file should be [SourceType::SourceName::file_name]:\n"
        data_export += "#1. Source Type - which can be library or history\n"
        data_export += "#2. Source Name - the name of the library or history.\n"
        data_export += "#3. Filename - specify the name of the file as it exists in the library or history.\n\n\n"
        data_export += "########################################\n\n\n"
        data_export += "### METADATA\n"
        data_export += "#######################################\n"
        data_export += "Workflow Name\t%s\n" % (stored.name)
        data_export += "Workflow id\t%s\n" % (trans.security.encode_id(workflow.id))
        data_export += "Project Name\t%s\n" % ("<Your_project_name>")
        data_export += "#######################################\n\n\n"

        data_export += "###TABLE DATA\n"
        data_export += "#######################################\n"
        data_export += "SampleName\t"

#        for attr in dir(stored):
#            print "obj.%s = %s" % (attr, getattr(stored, attr))

        #try: # use a try/finally block to restore the user's current history
            # Prepare each step
        missing_tools = []

        for step in workflow.steps:
                step_annotation = self.get_item_annotation_obj( trans.sa_session, trans.user, step )
                annotation_str = ""
                if step_annotation:
                    annotation_str = step_annotation.annotation

                #for property, value in vars(step).iteritems():
                #    print property, ": ", value

                step_input_connections = []
                for input_connection in step.input_connections:
                    step_input_connections.append(input_connection.input_name)
                    #print "BYE: %s" % input_connection
                #print "TOOL: %s" % step.tool_id
                #print "ALLO: %s" % step_input_connections

                step.upgrade_messages = {}

                # Contruct modules
                if step.type == 'tool' or step.type is None:
                    # Restore the tool state for the step
                    step.module = module_factory.from_workflow_step( trans, step )

                    # tool object
                    param_types = {}
                    if step.type == 'tool':
                        tool = trans.app.toolbox.get_tool( step.tool_id )
                        #print "TOOL INPUTS OBJ: %s" % tool.inputs

                        for param in tool.input_params:
                            #print param.name
                            if type(param) == DataToolParameter:
                                param_types[param.name] = type(param)
                        #print param_types

                    for key, value in step.module.state.inputs.items():
                        #print "KEY: %s: %s" % (key, value)
                        input_type = type(value)
                        #print input_type
                        ext = "%s" % (key)
                        if input_type == RuntimeValue and ext not in step_input_connections and key not in param_types:
                            #print "EXT: %s - COMP: %s" % (ext, step_input_connections)
                            data_export += annotation_str + "##Param::" +  str(step.order_index) + "::" + step.tool_id + "::" +  key + "\t"
                        elif input_type == dict:
#                            print "MAIN DICT"
#                            print value
                            data_export += _check_subdict(key, value, step, None, annotation_str, step_input_connections, param_types)
                        elif input_type == list:
                            for item in value:
                                #print "ITEM: %s" % item
                                if '__index__' in item:
                                    key_name = "%s_%s" % (key, item['__index__'])
                                else:
                                    key_name = "%s" % (key)
                                #print "KEY_NAME: %s" % key_name
                                data_export += _check_subdict(key_name, item, step, None, annotation_str, step_input_connections, param_types)
                            #print "LIST"
                            #print value
                            #for item in value:
                            #    print item


                    if not step.module:
                        if step.tool_id not in missing_tools:
                            missing_tools.append(step.tool_id)
                        continue
                    step.upgrade_messages = step.module.check_and_update_state()
                    if step.upgrade_messages:
                        has_upgrade_messages = True
                    # Any connected input needs to have value DummyDataset (these
                    # are not persisted so we need to do it every time)
                    step.module.add_dummy_datasets( connections=step.input_connections )
                    # Store state with the step
                    step.state = step.module.state
                    # Error dict
                    if step.tool_errors:
                        has_errors = True
                        errors[step.id] = step.tool_errors
                else:
                    ## Non-tool specific stuff?
                    step.module = module_factory.from_workflow_step( trans, step )
                    step.state = step.module.get_runtime_state()
                    #print step.__dict__['tool_inputs']['name']

                    #for attr in dir(step):
                    #    print "obj.%s = %s" % (attr, getattr(step, attr))

                    #print "%s: %s" % (step.module.name, step.__dict__['tool_inputs']['name'])
                    #data_export += annotation_str + "##SourceType::SourceName::%s\t" % (step.__dict__['tool_inputs']['name'])
                    data_export += annotation_str + "##SourceType::SourceName::%s\t" % (step.__dict__['label'])

        if missing_tools:
            stored.annotation = self.get_item_annotation_str( trans.sa_session, trans.user, stored )
            return trans.fill_template("workflow/run.mako",
                                       steps=[],
                                       workflow=stored,
                                       hide_fixed_params=hide_fixed_params,
                                       missing_tools = missing_tools)

        sname = stored.name
        valid_chars = '.,^_-()[]0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
        sname = ''.join(c in valid_chars and c or '_' for c in sname)[0:150]
        trans.response.headers["Content-Disposition"] = 'attachment; filename="Galaxy-API-Workflow-%s.txt"' % ( sname )
        trans.response.set_content_type( 'text/plain' )
        #print "DATA: %s" % data_export
        return data_export

def _check_subdict(key, value, step, existing_key, annotation_str, step_input_connections, param_types):
    data = ""
    if not isinstance(value, basestring) and value.__class__ != LibraryDatasetDatasetAssociation:
      #print "SUB_DICT: %s :: %s" %  (key, value)
      for partname, partval in value.items():
        if type ( partval ) == RuntimeValue:
          if existing_key is not None:
              ext = "%s|%s" % (existing_key, partname)
          else:
              ext = "%s|%s" % (key, partname)
          #print "EXT: %s - COMP: %s" % (ext, step_input_connections)
          #print "PARAM_TYPE: %s" % param_types
          if ext not in step_input_connections and partname not in param_types:
            ext = ext.replace("|", "::")
            #print "HELLO %s : %s" % (step.tool_id, partname)
            #print "LINE: " + annotation_str + "##Param::" + str(step.order_index) + "::" + step.tool_id + "::" +  ext + "\t"
#            data += annotation_str + "##Param::" + str(step.order_index) + "::" + step.tool_id + "::" +  key + "::" + partname + "\t"
            data += annotation_str + "##Param::" + str(step.order_index) + "::" + step.tool_id + "::" +  ext + "\t"
        elif isinstance(partval, types.NoneType):
            if existing_key:
                key_name = "%s|%s" % (existing_key, partname)
            else:
                key_name = "%s|%s" % (key, partname)
        elif type ( partval ) == dict:
            #print "DICT"
            #print partval
            key_name = "%s|%s" % (key, partname)
            #print key_name
            data += _check_subdict(partname, partval, step, key_name, annotation_str, step_input_connections, param_types)
        elif type( partval ) == list:
            #print "LIST-SUB"
            #print partval
            for item in partval:
                #print "ITEM: %s" % item
                if item.__class__ != LibraryDatasetDatasetAssociation and '__index__' in item:
                    key_name = "%s|%s_%s" % (key, partname, item['__index__'])
                else:
                    key_name = "%s|%s" % (key, partname)
#                print key_name
                data += _check_subdict(partname, item, step, key_name, annotation_str, step_input_connections, param_types)
    return data

def _build_workflow_on_str(instance_ds_names):
    # Returns suffix for new histories based on multi input iteration
    num_multi_inputs = len(instance_ds_names)
    if num_multi_inputs == 0:
        return ""
    elif num_multi_inputs == 1:
        return " on %s" % instance_ds_names[0]
    else:
        return " on %s and %s" % (", ".join(instance_ds_names[0:-1]), instance_ds_names[-1])


def _expand_multiple_inputs(kwargs):
    (single_inputs, matched_multi_inputs, multiplied_multi_inputs) = _split_inputs(kwargs)

    # Build up every combination of inputs to be run together.
    input_combos = _extend_with_matched_combos(single_inputs, matched_multi_inputs)
    input_combos = _extend_with_multiplied_combos(input_combos, multiplied_multi_inputs)

    # Input name that are multiply specified
    multi_input_keys = list(matched_multi_inputs.keys()) + list(multiplied_multi_inputs.keys())

    for input_combo in input_combos:
        for key, value in input_combo.items():
            kwargs[key] = value
        yield (kwargs, multi_input_keys)


def _extend_with_matched_combos(single_inputs, multi_inputs):
    if len(multi_inputs) == 0:
        return [single_inputs]

    matched_multi_inputs = []

    first_multi_input_key = next(iter(multi_inputs.keys()))
    first_multi_value = multi_inputs.get(first_multi_input_key)

    for value in first_multi_value:
        new_inputs = _copy_and_extend_inputs(single_inputs, first_multi_input_key, value)
        matched_multi_inputs.append(new_inputs)

    for multi_input_key, multi_input_values in multi_inputs.items():
        if multi_input_key == first_multi_input_key:
            continue
        if len(multi_input_values) != len(first_multi_value):
            raise Exception("Failed to match up multi-select inputs, must select equal number of data files in each multiselect")
        for index, value in enumerate(multi_input_values):
            matched_multi_inputs[index][multi_input_key] = value
    return matched_multi_inputs


def _extend_with_multiplied_combos(input_combos, multi_inputs):
    combos = input_combos

    for multi_input_key, multi_input_value in multi_inputs.items():
        iter_combos = []

        for combo in combos:
            for input_value in multi_input_value:
                iter_combos.append(_copy_and_extend_inputs(combo, multi_input_key, input_value))

        combos = iter_combos
    return combos


def _copy_and_extend_inputs(inputs, key, value):
    new_inputs = dict(inputs)
    new_inputs[key] = value
    return new_inputs


def _split_inputs(kwargs):
    """
    """
    input_keys = [a for a in kwargs if a.endswith('|input')]
    single_inputs = {}
    matched_multi_inputs = {}
    multiplied_multi_inputs = {}
    for input_key in input_keys:
        input_val = kwargs[input_key]
        if isinstance(input_val, list):
            input_base = input_key[:-len("|input")]
            mode_key = "%s|multi_mode" % input_base
            mode = kwargs.get(mode_key, "matched")
            if mode == "matched":
                matched_multi_inputs[input_key] = input_val
            else:
                multiplied_multi_inputs[input_key] = input_val
        else:
            single_inputs[input_key] = input_val
    return (single_inputs, matched_multi_inputs, multiplied_multi_inputs)
