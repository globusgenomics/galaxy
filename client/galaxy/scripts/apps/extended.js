/* global Galaxy, $ */
import "bootstrap";
export { GalaxyApp } from "galaxy";
import WorkflowView from "mvc/workflow/workflow-view";
import { Trackster } from "viz/trackster";
export { default as Trackster } from "viz/trackster";
import Circster from "viz/circster";
export { PhylovizView as phyloviz } from "viz/phyloviz";
export { SweepsterVisualization, SweepsterVisualizationView } from "viz/sweepster";
import GalaxyLibrary from "galaxy.library";
import AdminToolshed from "admin.toolshed";
import Masthead from "layout/masthead";
import user from "mvc/user/user-model";
import Modal from "mvc/ui/ui-modal";
export { default as pages } from "galaxy.pages";
export { createTabularDatasetChunkedView } from "mvc/dataset/data";
export { History } from "mvc/history/history-model";
export { HistoryContents } from "mvc/history/history-contents";
import MultiPanel from "mvc/history/multi-panel";
export { historyEntry as history } from "mvc/history/history-view";
export { default as HistoryViewAnnotated } from "mvc/history/history-view-annotated";
export { default as HistoryCopyDialog } from "mvc/history/copy-dialog";
export { default as HDAListItemEdit } from "mvc/history/hda-li-edit";
export { default as HDAModel } from "mvc/history/hda-model";
import addLogging from "utils/add-logging";
export { default as LegacyGridView } from "legacy/grid/grid-view";
export { default as run_stats } from "reports/run_stats";
export { default as ToolshedGroups } from "toolshed/toolshed.groups";

export { chart, chartUtilities } from "./chart";

// TODO: update this when galaxy singleton code is merged
if (window.Galaxy && window.Galaxy.debug === undefined) {
    //TODO: (kind of a temporary hack?) Must have Galaxy.logging for some of the imports
    //here; remove when imports are all fixed.
    addLogging(window.Galaxy, "GalaxyApp");
}

export function masthead(options) {
    if (!Galaxy.user) {
        Galaxy.user = new user.User(options.user_json);
    }
    if (!Galaxy.masthead) {
        Galaxy.masthead = new Masthead.View(options);
        Galaxy.modal = new Modal.View();
        $("#masthead").replaceWith(Galaxy.masthead.render().$el);
    }
}

export function adminToolshed(options) {
    new AdminToolshed.GalaxyApp(options);
}

export function trackster(options) {
    new Trackster.GalaxyApp(options);
}

export function circster(options) {
    new Circster.GalaxyApp(options);
}

export function workflow(options) {
    new WorkflowView(options);
}

export function library(options) {
    new GalaxyLibrary.GalaxyApp(options);
}

export function multiHistory(options) {
    let histories = new History.HistoryCollection([], {
        includeDeleted: options.includingDeleted,
        order: options.order,
        limitOnFirstFetch: options.limit,
        limitPerFetch: options.limit,
        currentHistoryId: options.current_history_id
    });
    let multipanel = new MultiPanel.MultiPanelColumns({
        el: $("#center").get(0),
        histories: histories
    });

    histories.fetchFirst({ silent: true }).done(function() {
        multipanel.createColumns();
        multipanel.render(0);
    });
}
