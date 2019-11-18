##
## Base template for sharing or publishing an item. Template expects the following parameters:
## (a) item - item to be shared.
##
<%!
    def inherit(context):
        if context.get('use_panels', False) == True:
            if context.get('webapp'):
                webapp = context.get('webapp')
            else:
                webapp = 'galaxy'
            return '/webapps/%s/base_panels.mako' % webapp
        else:
            return '/base.mako'
%>
<%inherit file="${inherit(context)}"/>

<%namespace file="/display_common.mako" import="*" />
<%namespace file="/message.mako" import="render_msg" />

##
## Page methods.
##

<%def name="init()">
<%
    self.has_left_panel=False
    self.has_right_panel=False
    self.message_box_visible=False
    self.overlay_visible=False
    self.message_box_class=""
    #self.active_view=""
    self.body_class=""
%>
</%def>

<%def name="title()">
    Exporting workflow parameters for batch submission ${get_class_display_name( item.__class__ )} '${get_item_name( item )}'
</%def>

<%def name="javascripts()">
    ${parent.javascripts()}
    <script type="text/javascript">
    $(document).ready( function() {
        //
        // Set up slug-editing functionality.
        //
        var on_start = function( text_elt ) {
            // Replace URL with URL text.
            $('#item-url').hide();
            $('#item-url-text').show();
            
            // Allow only lowercase alphanumeric and '-' characters in slug.
            text_elt.keyup(function(){
                text_elt.val( $(this).val().replace(/\s+/g,'-').replace(/[^a-zA-Z0-9\-]/g,'').toLowerCase() )
            });
        };
        
        var on_finish = function( text_elt ) {
            // Replace URL text with URL.
            $('#item-url-text').hide();
            $('#item-url').show();
            
            // Set URL to new value.
            var new_url = $('#item-url-text').text();
            var item_url_obj = $('#item-url');
            item_url_obj.attr( "href", new_url );
            item_url_obj.text( new_url );
        };
        
        <% controller_name = get_controller_name( item ) %>
        async_save_text("edit-identifier", "item-identifier", "${h.url_for( controller=controller_name, action='set_slug_async', id=trans.security.encode_id( item.id ) )}", "new_slug", null, false, 0, on_start, on_finish); 
    });
    </script>
</%def>

<%def name="stylesheets()">
    ${parent.stylesheets()}
    <style>
        ## Put some whitespace before each section header.
        h3 {
            margin-top: 2em;
        }
        input.action-button {
            margin-left: 0;
        }
        ## If page is displayed in panels, pad from edges for readability.
        %if context.get('use_panels'):
            div#center {
	        overflow: auto;
                padding: 10px;
            }
        %endif
    </style>
</%def>

<%def name="center_panel()">
    ${self.body()}
</%def>

<%def name="body()">
    ## Set use_panels var for use in page's URLs.
    <% use_panels = context.get('use_panels', False)  %>

    ## Render message.
    %if message:
        ${render_msg( message, status )}
    %endif

    <%
        #
        # Setup and variables needed for page.
        #
    
        # Get class name strings.
        item_class_name = get_class_display_name( item.__class__ ) 
        item_class_name_lc = item_class_name.lower()
        item_class_plural_name = get_class_plural_display_name( item.__class__ )
        item_class_plural_name_lc = item_class_plural_name.lower()
        
        # Get item name.
        item_name = get_item_name(item)
    %>

    <h2>Export Workflow Parameters for API Batch Submission:  ${item_class_name} '${item_name}'</h2>

    ## Require that user have a public username before anything.
    %if trans.get_user().username is None or trans.get_user().username is "":
        <p>To export parameters of a ${item_class_name_lc}, you must create a public username:</p>
        
        <form action="${h.url_for( action='set_public_username', id=trans.security.encode_id( item.id ) )}"     
                method="POST">
            <div class="form-row">
                <label>Public Username:</label>
                <div class="form-row-input">
                    <input type="text" name="username" size="40"/>
                </div>
            </div>
            <div style="clear: both"></div>
            <div class="form-row">
                <input class="action-button" type="submit" name="Set Username" value="Set Username"/>
            </div>
        </form>
    %else:
        ## User has a public username, so you can export parameters.
        <h3>Export parameters of ${item_class_name_lc} for API batch submission</h3>
    
            <div>
                <p>The Globus Genomics Galaxy API allows submission of a user defined workflow multiple times. You will need to do the following:</p>
		<ol>
		<li><b>API Key</b> - You will need to generate an API key to identify yourself with the Galaxy server. This will be a required input for the python script that will submit your jobs. If you don't have an API key, please generate one by following the <a href="http://wiki.galaxyproject.org/Learn/API#Enabling_the_API" target="_new">instructions</a>.</li> <br />
		<li><b>Workflow parameters table</b> - You can create a workflow through the workflow generator. You can set which parameters should be set at run time. You can download a skeleton of the table structure you will need to fill out with your input files and parameters that are specific to your workflow. Please don't modify any of the header rows or non-commented lines. You can download the file by clicking on the link below.</li>
                <a href="${h.url_for( controller=get_controller_name( item ), action='display_by_username_and_slug', username=item.user.username, slug=item.slug, format='export-parameters' )}">Export ${item_class_name} Parameters for batch submission</a><br />
                <div class="toolParamHelp">
	        <h3>Instructions</h3>
		<p>The table you have downloaded will help you submit a workflow in batch mode. A workflow can be as simple as running one tool (one process), to having hundreds of tools (hundreds of processes). The structure of the table file downloaded is shown in a figure below. The contents of the file can be separated into three sections:</p>
		<ul>
		<li>General instructions</li>
		<li>Workflow metadata - The only parameter the user should modify here is the "project name". All other items in the metadata section should be left alone or the submission process will not work</li>
		<li>Table data - the header will be supplied in the file. It's important to note that the header row cannot be modified by the user. If so, this will cause errors when submitting the workflow. After the header row, the user should feel free to input the necessary file location, names and parameters for their workflow. The items in the header that will be allowed to be modified are for:
		<ul>
		<li>Input files required for a tool</li>
		<li>Parameter values for tools which have been previously set by the user to be filled out at runtime</li>
		</ul>
		</li>
		</ul>
		<br /><br />
                <p>The following example should explain what is needed to run workflows in batch mode.</p>
		<ol>
		<li>Create an API key if you don't have one yet by following instructions <a href="http://wiki.galaxyproject.org/Learn/API#Enabling_the_API" target="_new">here</a>.
		<ul>
		<li>My API key on my local instance is: 3ab0cce215d049b1a1a2bdfa1eaf6b6b</li>
		<li>NOTE: This API key will not work on your instance, thus you will need to generate one so that it identifies you and creates the histories under your user account.</li>
		</ul>
		</li><br />
		<li>Create a workflow you wish to submit in batch mode. 
		<ul>
		<li><a href="http://dev.globusgenomics.org/u/arodri7/w/api-batch-test-worfklow" target="_new"">My example is a simple workflow</a> that takes two VCF files, cuts a few columns, adds a column with a user specified text, then concatenates the two VCF files.</li>
		<li>This workflow in available in the "Shared Data" and "Published Workflows" links in the main menu or <a href="http://dev.globusgenomics.org/u/arodri7/w/api-batch-test-worfklow" target="_new">here</a>. You should be able to import it to your user account by clicking on the "Import workflow" link.</li>
                <li>For each tool's input file you will need to create an "Input Dataset" box by searching for the label "Inputs" in the tool panel and selecting "Input Datasets". If this is not specified, you will not be able to specify an input file for that tool.</li>
		<li>Notice that the workflow has two inputs labeled as Ref1 and Ref2. Adding the "Input Dataset" box tool allowed this to happen.</li><image src="../static/images/API-test-workflow.png " />
		<li>In addition, the "add column" tool has the "Add this value" parameter which has been set to be filled out at runtime. This allows the user to see the column in the batch mode table you will soon download.</li><image src="../static/images/API-test-workflow-runtime.png " />
		</ul>
	    </li><br />
			<li>On the "Workflows" page click on the workflow you wish to submit in batch mode. A sub-menu will appear. Click on the "Submit via API batch mode" link. The link will take you to a page showing the same instructions you are following.</li><image src="../static/images/API-test-workflow-menu.png " />
                        <li> You will need to download a couple files:
			<ul>
			<li><a href="https://bitbucket.org/alex_rodriguez/api-workflow-batch-submission" target="_new">General workflow python submit script</a> - Please save this file on your computer and don't modify it. This script will allow you to interact with the BioBlend library to submit all the jobs to the Globus Genomics instance.</li>
			<li><a href="#param_table">Export Workflow Parameters for Batch Submission</a> - Clicking the button at the top of the page will download a tab-delimited file which can be opened with your favorite text editor.</li>
			</ul>
			</li><br />
			<li>Click the <a href="#param_table">"Export Workflow Parameters for Batch Submission"</a> button to download the file you will need to edit.</li><br />
			<li>Open the tab-delimited file you have downloaded. It should look similar to the figure below.</li><image src="../static/images/API-test-workflow-tabfile.png " /><br />
			<li>Modify the file similar to the figure above by adding the new rows shown. Each row represents a different submission of the workflow. Each submission will be run in parallel. Make sure the columns are separated by a tab character.</li><image src="../static/images/API-test-workflow-data.png " /><br />
			<li>Save the file as a text file. Any other type of file will not work.</li><br />
                        <li>Upload the file back to your Globus Genomics instance.</li><br />
                        <li>Towards the bottom of the tool panel, click on the "Batch Submit" tool and select the Table file for the workflow you want to submit.</li><br />
                        <li>Click on "Execute".</li><br />
                        <li>If you go to the instance and look at your saved histories, you will see one new history for each row you are submitting. The histories are named using a timestamp to make sure they have unique names.</li><image src="../static/images/API-test-workflow-history.png " />
        		</div>
                        <br />


    %endif

    <br /><br />
</%def>
