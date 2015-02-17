define(["utils/levenshtein","utils/natural-sort","mvc/base-mvc","utils/localization"],function(h,b,f,c){var i=Backbone.View.extend(f.LoggableMixin).extend({tagName:"li",className:"dataset paired",initialize:function(l){this.pair=l.pair||{}},template:_.template(['<span class="forward-dataset-name flex-column"><%= pair.forward.name %></span>','<span class="pair-name-column flex-column">','<span class="pair-name"><%= pair.name %></span>',"</span>",'<span class="reverse-dataset-name flex-column"><%= pair.reverse.name %></span>'].join("")),render:function(){this.$el.attr("draggable",true).data("pair",this.pair).html(this.template({pair:this.pair})).addClass("flex-column-container");return this},events:{dragstart:"_dragstart",dragend:"_dragend",dragover:"_sendToParent",drop:"_sendToParent"},_dragstart:function(l){l.currentTarget.style.opacity="0.4";if(l.originalEvent){l=l.originalEvent}l.dataTransfer.effectAllowed="move";l.dataTransfer.setData("text/plain",JSON.stringify(this.pair));this.$el.parent().trigger("pair.dragstart",[this])},_dragend:function(l){l.currentTarget.style.opacity="1.0";this.$el.parent().trigger("pair.dragend",[this])},_sendToParent:function(l){this.$el.parent().trigger(l)},toString:function(){return"PairView("+this.pair.name+")"}});function g(m){m=m||{};m.createPair=m.createPair||function l(r){this.debug("creating pair:",r.listA[r.indexA].name,r.listB[r.indexB].name);r=r||{};return this._pair(r.listA.splice(r.indexA,1)[0],r.listB.splice(r.indexB,1)[0],{silent:true})};var o=[];function q(){if(!o.length){o=[new RegExp(this.filters[0]),new RegExp(this.filters[1])]}return o}m.preprocessMatch=m.preprocessMatch||function n(s){var r=q.call(this);return _.extend(s,{matchTo:s.matchTo.name.replace(r[0],""),possible:s.possible.name.replace(r[1],"")})};return function p(t){this.debug("autopair _strategy ---------------------------");t=t||{};var r=t.listA,A=t.listB,z=0,y,v={score:0,index:null},x=[];this.debug("starting list lens:",r.length,A.length);this.debug("bestMatch (starting):",JSON.stringify(v,null,"  "));while(z<r.length){var w=r[z];v.score=0;for(y=0;y<A.length;y++){var u=A[y];this.debug(z+":"+w.name);this.debug(y+":"+u.name);if(r[z]!==A[y]){v=m.match.call(this,m.preprocessMatch.call(this,{matchTo:w,possible:u,index:y,bestMatch:v}));this.debug("bestMatch:",JSON.stringify(v,null,"  "));if(v.score===1){this.debug("breaking early due to perfect match");break}}}var s=m.scoreThreshold.call(this);this.debug("scoreThreshold:",s);this.debug("bestMatch.score:",v.score);if(v.score>=s){this.debug("creating pair");x.push(m.createPair.call(this,{listA:r,indexA:z,listB:A,indexB:v.index}));this.debug("list lens now:",r.length,A.length)}else{z+=1}if(!r.length||!A.length){return x}}this.debug("paired:",JSON.stringify(x,null,"  "));this.debug("autopair _strategy ---------------------------");return x}}var k=Backbone.View.extend(f.LoggableMixin).extend({className:"collection-creator flex-row-container",initialize:function(l){this.metric("PairedCollectionCreator.initialize",l);l=_.defaults(l,{datasets:[],filters:this.DEFAULT_FILTERS,automaticallyPair:true,strategy:"lcs",matchPercentage:0.9,twoPassAutopairing:true});this.initialList=l.datasets;this.historyId=l.historyId;this.filters=this.commonFilters[l.filters]||this.commonFilters[this.DEFAULT_FILTERS];if(_.isArray(l.filters)){this.filters=l.filters}this.automaticallyPair=l.automaticallyPair;this.strategy=this.strategies[l.strategy]||this.strategies[this.DEFAULT_STRATEGY];if(_.isFunction(l.strategy)){this.strategy=l.strategy}this.matchPercentage=l.matchPercentage;this.twoPassAutopairing=l.twoPassAutopairing;this.removeExtensions=true;this.oncancel=l.oncancel;this.oncreate=l.oncreate;this.autoscrollDist=l.autoscrollDist||24;this.unpairedPanelHidden=false;this.pairedPanelHidden=false;this.$dragging=null;this._setUpBehaviors();this._dataSetUp()},commonFilters:{illumina:["_1","_2"],Rs:["_R1","_R2"]},DEFAULT_FILTERS:"illumina",strategies:{simple:"autopairSimple",lcs:"autopairLCS",levenshtein:"autopairLevenshtein"},DEFAULT_STRATEGY:"lcs",_dataSetUp:function(){this.paired=[];this.unpaired=[];this.selectedIds=[];this._sortInitialList();this._ensureIds();this.unpaired=this.initialList.slice(0);if(this.automaticallyPair){this.autoPair();this.once("rendered:initial",function(){this.trigger("autopair")})}},_sortInitialList:function(){this._sortDatasetList(this.initialList)},_sortDatasetList:function(l){l.sort(function(n,m){return b(n.name,m.name)});return l},_ensureIds:function(){this.initialList.forEach(function(l){if(!l.hasOwnProperty("id")){l.id=_.uniqueId()}});return this.initialList},_splitByFilters:function(){var o=this.filters.map(function(p){return new RegExp(p)}),m=[[],[]];function n(p,q){return q.test(p.name)}this.unpaired.forEach(function l(p){o.forEach(function(r,q){if(n(p,r)){m[q].push(p)}})});return m},_addToUnpaired:function(m){var l=function(n,p){if(n===p){return n}var o=Math.floor((p-n)/2)+n,q=b(m.name,this.unpaired[o].name);if(q<0){return l(n,o)}else{if(q>0){return l(o+1,p)}}while(this.unpaired[o]&&this.unpaired[o].name===m.name){o++}return o}.bind(this);this.unpaired.splice(l(0,this.unpaired.length),0,m)},autoPair:function(n){var m=this._splitByFilters(),l=[];if(this.twoPassAutopairing){l=this.autopairSimple({listA:m[0],listB:m[1]});m=this._splitByFilters()}n=n||this.strategy;m=this._splitByFilters();l=l.concat(this[n].call(this,{listA:m[0],listB:m[1]}));return l},autopairSimple:g({scoreThreshold:function(){return 1},match:function j(l){l=l||{};if(l.matchTo===l.possible){return{index:l.index,score:1}}return l.bestMatch}}),autopairLevenshtein:g({scoreThreshold:function(){return this.matchPercentage},match:function e(m){m=m||{};var n=h(m.matchTo,m.possible),l=1-(n/(Math.max(m.matchTo.length,m.possible.length)));if(l>m.bestMatch.score){return{index:m.index,score:l}}return m.bestMatch}}),autopairLCS:g({scoreThreshold:function(){return this.matchPercentage},match:function e(n){n=n||{};var l=this._naiveStartingAndEndingLCS(n.matchTo,n.possible).length,m=l/(Math.max(n.matchTo.length,n.possible.length));if(m>n.bestMatch.score){return{index:n.index,score:m}}return n.bestMatch}}),_naiveStartingAndEndingLCS:function(o,m){var p="",q="",n=0,l=0;while(n<o.length&&n<m.length){if(o[n]!==m[n]){break}p+=o[n];n+=1}if(n===o.length){return o}if(n===m.length){return m}n=(o.length-1);l=(m.length-1);while(n>=0&&l>=0){if(o[n]!==m[l]){break}q=[o[n],q].join("");n-=1;l-=1}return p+q},_pair:function(n,l,m){m=m||{};var o=this._createPair(n,l,m.name);this.paired.push(o);this.unpaired=_.without(this.unpaired,n,l);if(!m.silent){this.trigger("pair:new",o)}return o},_createPair:function(n,l,m){if(!(n&&l)||(n===l)){throw new Error("Bad pairing: "+[JSON.stringify(n),JSON.stringify(l)])}m=m||this._guessNameForPair(n,l);return{forward:n,name:m,reverse:l}},_guessNameForPair:function(o,m,q){q=(q!==undefined)?(q):(this.removeExtensions);var l=o.name,p=m.name,n=this._naiveStartingAndEndingLCS(l.replace(this.filters[0],""),p.replace(this.filters[1],""));if(q){var r=n.lastIndexOf(".");if(r>0){var s=n.slice(r,n.length);n=n.replace(s,"");l=l.replace(s,"");p=p.replace(s,"")}}return n||(l+" & "+p)},_unpair:function(m,l){l=l||{};if(!m){throw new Error("Bad pair: "+JSON.stringify(m))}this.paired=_.without(this.paired,m);this._addToUnpaired(m.forward);this._addToUnpaired(m.reverse);if(!l.silent){this.trigger("pair:unpair",[m])}return m},unpairAll:function(){var l=[];while(this.paired.length){l.push(this._unpair(this.paired[0],{silent:true}))}this.trigger("pair:unpair",l)},_pairToJSON:function(m,l){l=l||"hda";return{collection_type:"paired",src:"new_collection",name:m.name,element_identifiers:[{name:"forward",id:m.forward.id,src:l},{name:"reverse",id:m.reverse.id,src:l}]}},createList:function(n){var o=this,m="/api/histories/"+this.historyId+"/contents/dataset_collections";var l={type:"dataset_collection",collection_type:"list:paired",name:_.escape(n||o.$(".collection-name").val()),element_identifiers:o.paired.map(function(p){return o._pairToJSON(p)})};return jQuery.ajax(m,{type:"POST",contentType:"application/json",dataType:"json",data:JSON.stringify(l)}).fail(function(r,p,q){o._ajaxErrHandler(r,p,q)}).done(function(p,q,r){o.trigger("collection:created",p,q,r);o.metric("collection:created",p);if(typeof o.oncreate==="function"){o.oncreate.call(this,p,q,r)}})},_ajaxErrHandler:function(o,l,n){this.error(o,l,n);var m=c("An error occurred while creating this collection");if(o){if(o.readyState===0&&o.status===0){m+=": "+c("Galaxy could not be reached and may be updating.")+c(" Try again in a few minutes.")}else{if(o.responseJSON){m+="<br /><pre>"+JSON.stringify(o.responseJSON)+"</pre>"}else{m+=": "+n}}}creator._showAlert(m,"alert-danger")},render:function(l,m){this.$el.empty().html(k.templates.main());this._renderHeader(l);this._renderMiddle(l);this._renderFooter(l);this._addPluginComponents();this.trigger("rendered",this);return this},_renderHeader:function(m,n){var l=this.$(".header").empty().html(k.templates.header()).find(".help-content").prepend($(k.templates.helpContent()));this._renderFilters();return l},_renderFilters:function(){return this.$(".forward-column .column-header input").val(this.filters[0]).add(this.$(".reverse-column .column-header input").val(this.filters[1]))},_renderMiddle:function(m,n){var l=this.$(".middle").empty().html(k.templates.middle());if(this.unpairedPanelHidden){this.$(".unpaired-columns").hide()}else{if(this.pairedPanelHidden){this.$(".paired-columns").hide()}}this._renderUnpaired();this._renderPaired();return l},_renderUnpaired:function(q,r){var o=this,p,m,l=[],n=this._splitByFilters();this.$(".forward-column .title").text([n[0].length,c("unpaired forward")].join(" "));this.$(".forward-column .unpaired-info").text(this._renderUnpairedDisplayStr(this.unpaired.length-n[0].length));this.$(".reverse-column .title").text([n[1].length,c("unpaired reverse")].join(" "));this.$(".reverse-column .unpaired-info").text(this._renderUnpairedDisplayStr(this.unpaired.length-n[1].length));this.$(".unpaired-columns .column-datasets").empty();this.$(".autopair-link").toggle(this.unpaired.length!==0);if(this.unpaired.length===0){this._renderUnpairedEmpty();return}m=n[1].map(function(t,s){if((n[0][s]!==undefined)&&(n[0][s]!==t)){l.push(o._renderPairButton())}return o._renderUnpairedDataset(t)});p=n[0].map(function(s){return o._renderUnpairedDataset(s)});if(!p.length&&!m.length){this._renderUnpairedNotShown();return}this.$(".unpaired-columns .forward-column .column-datasets").append(p).add(this.$(".unpaired-columns .paired-column .column-datasets").append(l)).add(this.$(".unpaired-columns .reverse-column .column-datasets").append(m));this._adjUnpairedOnScrollbar()},_renderUnpairedDisplayStr:function(l){return["(",l," ",c("filtered out"),")"].join("")},_renderUnpairedDataset:function(l){return $("<li/>").attr("id","dataset-"+l.id).addClass("dataset unpaired").attr("draggable",true).addClass(l.selected?"selected":"").append($("<span/>").addClass("dataset-name").text(l.name)).data("dataset",l)},_renderPairButton:function(){return $("<li/>").addClass("dataset unpaired").append($("<span/>").addClass("dataset-name").text(c("Pair these datasets")))},_renderUnpairedEmpty:function(){var l=$('<div class="empty-message"></div>').text("("+c("no remaining unpaired datasets")+")");this.$(".unpaired-columns .paired-column .column-datasets").empty().prepend(l);return l},_renderUnpairedNotShown:function(){var l=$('<div class="empty-message"></div>').text("("+c("no datasets were found matching the current filters")+")");this.$(".unpaired-columns .paired-column .column-datasets").empty().prepend(l);return l},_adjUnpairedOnScrollbar:function(){var o=this.$(".unpaired-columns").last(),p=this.$(".unpaired-columns .reverse-column .dataset").first();if(!p.size()){return}var l=o.offset().left+o.outerWidth(),n=p.offset().left+p.outerWidth(),m=Math.floor(l)-Math.floor(n);this.$(".unpaired-columns .forward-column").css("margin-left",(m>0)?m:0)},_renderPaired:function(m,n){this.$(".paired-column-title .title").text([this.paired.length,c("paired")].join(" "));this.$(".unpair-all-link").toggle(this.paired.length!==0);if(this.paired.length===0){this._renderPairedEmpty();return}else{this.$(".remove-extensions-link").show()}this.$(".paired-columns .column-datasets").empty();var l=this;this.paired.forEach(function(q,o){var p=new i({pair:q});l.$(".paired-columns .column-datasets").append(p.render().$el).append(['<button class="unpair-btn">','<span class="fa fa-unlink" title="',c("Unpair"),'"></span>',"</button>"].join(""))})},_renderPairedEmpty:function(){var l=$('<div class="empty-message"></div>').text("("+c("no paired datasets yet")+")");this.$(".paired-columns .column-datasets").empty().prepend(l);return l},_renderFooter:function(m,n){var l=this.$(".footer").empty().html(k.templates.footer());this.$(".remove-extensions").prop("checked",this.removeExtensions);if(typeof this.oncancel==="function"){this.$(".cancel-create.btn").show()}return l},_addPluginComponents:function(){this._chooseFiltersPopover(".choose-filters-link");this.$(".help-content i").hoverhighlight(".collection-creator","rgba( 64, 255, 255, 1.0 )")},_chooseFiltersPopover:function(l){function m(p,o){return['<button class="filter-choice btn" ','data-forward="',p,'" data-reverse="',o,'">',c("Forward"),": ",p,", ",c("Reverse"),": ",o,"</button>"].join("")}var n=$(_.template(['<div class="choose-filters">','<div class="help">',c("Choose from the following filters to change which unpaired reads are shown in the display"),":</div>",_.values(this.commonFilters).map(function(o){return m(o[0],o[1])}).join(""),"</div>"].join(""))({}));return this.$(l).popover({container:".collection-creator",placement:"bottom",html:true,content:n})},_validationWarning:function(m,l){var n="validation-warning";if(m==="name"){m=this.$(".collection-name").add(this.$(".collection-name-prompt"));this.$(".collection-name").focus().select()}if(l){m=m||this.$("."+n);m.removeClass(n)}else{m.addClass(n)}},_setUpBehaviors:function(){this.once("rendered",function(){this.trigger("rendered:initial",this)});this.on("pair:new",function(){this._renderUnpaired();this._renderPaired();this.$(".paired-columns").scrollTop(8000000)});this.on("pair:unpair",function(l){this._renderUnpaired();this._renderPaired();this.splitView()});this.on("filter-change",function(){this.filters=[this.$(".forward-unpaired-filter input").val(),this.$(".reverse-unpaired-filter input").val()];this.metric("filter-change",this.filters);this._renderFilters();this._renderUnpaired()});this.on("autopair",function(){this._renderUnpaired();this._renderPaired();var l,m=null;if(this.paired.length){m="alert-success";l=this.paired.length+" "+c("pairs created");if(!this.unpaired.length){l+=": "+c("all datasets have been successfully paired");this.hideUnpaired();this.$(".collection-name").focus()}}else{l=c("Could not automatically create any pairs from the given dataset names")}this._showAlert(l,m)});return this},events:{"click .more-help":"_clickMoreHelp","click .less-help":"_clickLessHelp","click .header .alert button":"_hideAlert","click .forward-column .column-title":"_clickShowOnlyUnpaired","click .reverse-column .column-title":"_clickShowOnlyUnpaired","click .unpair-all-link":"_clickUnpairAll","change .forward-unpaired-filter input":function(l){this.trigger("filter-change")},"focus .forward-unpaired-filter input":function(l){$(l.currentTarget).select()},"click .autopair-link":"_clickAutopair","click .choose-filters .filter-choice":"_clickFilterChoice","click .clear-filters-link":"_clearFilters","change .reverse-unpaired-filter input":function(l){this.trigger("filter-change")},"focus .reverse-unpaired-filter input":function(l){$(l.currentTarget).select()},"click .forward-column .dataset.unpaired":"_clickUnpairedDataset","click .reverse-column .dataset.unpaired":"_clickUnpairedDataset","click .paired-column .dataset.unpaired":"_clickPairRow","click .unpaired-columns":"clearSelectedUnpaired","mousedown .unpaired-columns .dataset":"_mousedownUnpaired","click .paired-column-title":"_clickShowOnlyPaired","mousedown .flexible-partition-drag":"_startPartitionDrag","click .paired-columns .dataset.paired":"selectPair","click .paired-columns":"clearSelectedPaired","click .paired-columns .pair-name":"_clickPairName","click .unpair-btn":"_clickUnpair","dragover .paired-columns .column-datasets":"_dragoverPairedColumns","drop .paired-columns .column-datasets":"_dropPairedColumns","pair.dragstart .paired-columns .column-datasets":"_pairDragstart","pair.dragend   .paired-columns .column-datasets":"_pairDragend","change .remove-extensions":function(l){this.toggleExtensions()},"change .collection-name":"_changeName","keydown .collection-name":"_nameCheckForEnter","click .cancel-create":function(l){if(typeof this.oncancel==="function"){this.oncancel.call(this)}},"click .create-collection":"_clickCreate"},_clickMoreHelp:function(l){this.$(".main-help").addClass("expanded");this.$(".more-help").hide()},_clickLessHelp:function(l){this.$(".main-help").removeClass("expanded");this.$(".more-help").show()},_showAlert:function(m,l){l=l||"alert-danger";this.$(".main-help").hide();this.$(".header .alert").attr("class","alert alert-dismissable").addClass(l).show().find(".alert-message").html(m)},_hideAlert:function(l){this.$(".main-help").show();this.$(".header .alert").hide()},_clickShowOnlyUnpaired:function(l){if(this.$(".paired-columns").is(":visible")){this.hidePaired()}else{this.splitView()}},_clickShowOnlyPaired:function(l){if(this.$(".unpaired-columns").is(":visible")){this.hideUnpaired()}else{this.splitView()}},hideUnpaired:function(l,m){this.unpairedPanelHidden=true;this.pairedPanelHidden=false;this._renderMiddle(l,m)},hidePaired:function(l,m){this.unpairedPanelHidden=false;this.pairedPanelHidden=true;this._renderMiddle(l,m)},splitView:function(l,m){this.unpairedPanelHidden=this.pairedPanelHidden=false;this._renderMiddle(l,m);return this},_clickUnpairAll:function(l){this.metric("unpairAll");this.unpairAll()},_clickAutopair:function(m){var l=this.autoPair();this.metric("autopair",l.length,this.unpaired.length);this.trigger("autopair")},_clickFilterChoice:function(m){var l=$(m.currentTarget);this.$(".forward-unpaired-filter input").val(l.data("forward"));this.$(".reverse-unpaired-filter input").val(l.data("reverse"));this._hideChooseFilters();this.trigger("filter-change")},_hideChooseFilters:function(){this.$(".choose-filters-link").popover("hide");this.$(".popover").css("display","none")},_clearFilters:function(l){this.$(".forward-unpaired-filter input").val("");this.$(".reverse-unpaired-filter input").val("");this.trigger("filter-change")},_clickUnpairedDataset:function(l){l.stopPropagation();return this.toggleSelectUnpaired($(l.currentTarget))},toggleSelectUnpaired:function(n,m){m=m||{};var o=n.data("dataset"),l=m.force!==undefined?m.force:!n.hasClass("selected");if(!n.size()||o===undefined){return n}if(l){n.addClass("selected");if(!m.waitToPair){this.pairAllSelected()}}else{n.removeClass("selected")}return n},pairAllSelected:function(m){m=m||{};var n=this,o=[],l=[],p=[];n.$(".unpaired-columns .forward-column .dataset.selected").each(function(){o.push($(this).data("dataset"))});n.$(".unpaired-columns .reverse-column .dataset.selected").each(function(){l.push($(this).data("dataset"))});o.length=l.length=Math.min(o.length,l.length);o.forEach(function(r,q){try{p.push(n._pair(r,l[q],{silent:true}))}catch(s){n.error(s)}});if(p.length&&!m.silent){this.trigger("pair:new",p)}return p},clearSelectedUnpaired:function(){this.$(".unpaired-columns .dataset.selected").removeClass("selected")},_mousedownUnpaired:function(n){if(n.shiftKey){var m=this,l=$(n.target).addClass("selected"),o=function(p){m.$(p.target).filter(".dataset").addClass("selected")};l.parent().on("mousemove",o);$(document).one("mouseup",function(p){l.parent().off("mousemove",o);m.pairAllSelected()})}},_clickPairRow:function(n){var o=$(n.currentTarget).index(),m=$(".unpaired-columns .forward-column .dataset").eq(o).data("dataset"),l=$(".unpaired-columns .reverse-column .dataset").eq(o).data("dataset");this._pair(m,l)},_startPartitionDrag:function(m){var l=this,p=m.pageY;$("body").css("cursor","ns-resize");l.$(".flexible-partition-drag").css("color","black");function o(q){l.$(".flexible-partition-drag").css("color","");$("body").css("cursor","").unbind("mousemove",n)}function n(q){var r=q.pageY-p;if(!l.adjPartition(r)){$("body").trigger("mouseup")}l._adjUnpairedOnScrollbar();p+=r}$("body").mousemove(n);$("body").one("mouseup",o)},adjPartition:function(m){var l=this.$(".unpaired-columns"),n=this.$(".paired-columns"),o=parseInt(l.css("height"),10),p=parseInt(n.css("height"),10);o=Math.max(10,o+m);p=p-m;var q=m<0;if(q){if(this.unpairedPanelHidden){return false}else{if(o<=10){this.hideUnpaired();return false}}}else{if(this.unpairedPanelHidden){l.show();this.unpairedPanelHidden=false}}if(!q){if(this.pairedPanelHidden){return false}else{if(p<=15){this.hidePaired();return false}}}else{if(this.pairedPanelHidden){n.show();this.pairedPanelHidden=false}}l.css({height:o+"px",flex:"0 0 auto"});return true},selectPair:function(l){l.stopPropagation();$(l.currentTarget).toggleClass("selected")},clearSelectedPaired:function(l){this.$(".paired-columns .dataset.selected").removeClass("selected")},_clickPairName:function(o){o.stopPropagation();var q=$(o.currentTarget),n=q.parent().parent(),m=n.index(".dataset.paired"),p=this.paired[m],l=prompt("Enter a new name for the pair:",p.name);if(l){p.name=l;p.customizedName=true;q.text(p.name)}},_clickUnpair:function(m){var l=Math.floor($(m.currentTarget).index(".unpair-btn"));this._unpair(this.paired[l])},_dragoverPairedColumns:function(o){o.preventDefault();var m=this.$(".paired-columns .column-datasets");this._checkForAutoscroll(m,o.originalEvent.clientY);var n=this._getNearestPairedDatasetLi(o.originalEvent.clientY);$(".paired-drop-placeholder").remove();var l=$('<div class="paired-drop-placeholder"></div>');if(!n.size()){m.append(l)}else{n.before(l)}},_checkForAutoscroll:function(l,r){var p=2;var q=l.offset(),o=l.scrollTop(),m=r-q.top,n=(q.top+l.outerHeight())-r;if(m>=0&&m<this.autoscrollDist){l.scrollTop(o-p)}else{if(n>=0&&n<this.autoscrollDist){l.scrollTop(o+p)}}},_getNearestPairedDatasetLi:function(r){var o=4,m=this.$(".paired-columns .column-datasets li").toArray();for(var n=0;n<m.length;n++){var q=$(m[n]),p=q.offset().top,l=Math.floor(q.outerHeight()/2)+o;if(p+l>r&&p-l<r){return q}}return $()},_dropPairedColumns:function(m){m.preventDefault();m.dataTransfer.dropEffect="move";var l=this._getNearestPairedDatasetLi(m.originalEvent.clientY);if(l.size()){this.$dragging.insertBefore(l)}else{this.$dragging.insertAfter(this.$(".paired-columns .unpair-btn").last())}this._syncPairsToDom();return false},_syncPairsToDom:function(){var l=[];this.$(".paired-columns .dataset.paired").each(function(){l.push($(this).data("pair"))});this.paired=l;this._renderPaired()},_pairDragstart:function(m,n){n.$el.addClass("selected");var l=this.$(".paired-columns .dataset.selected");this.$dragging=l},_pairDragend:function(l,m){$(".paired-drop-placeholder").remove();this.$dragging=null},toggleExtensions:function(m){var l=this;l.removeExtensions=(m!==undefined)?(m):(!l.removeExtensions);_.each(l.paired,function(n){if(n.customizedName){return}n.name=l._guessNameForPair(n.forward,n.reverse)});l._renderPaired();l._renderFooter()},_changeName:function(l){this._validationWarning("name",!!this._getName())},_nameCheckForEnter:function(l){if(l.keyCode===13){this._clickCreate()}},_getName:function(){return _.escape(this.$(".collection-name").val())},_clickCreate:function(m){var l=this._getName();if(!l){this._validationWarning("name")}else{this.createList()}},_printList:function(m){var l=this;_.each(m,function(n){if(m===l.paired){l._printPair(n)}else{}})},_printPair:function(l){this.debug(l.forward.name,l.reverse.name,": ->",l.name)},toString:function(){return"PairedCollectionCreator"}});k.templates=k.templates||{main:_.template(['<div class="header flex-row no-flex"></div>','<div class="middle flex-row flex-row-container"></div>','<div class="footer flex-row no-flex">'].join("")),header:_.template(['<div class="main-help well clear">','<a class="more-help" href="javascript:void(0);">',c("More help"),"</a>",'<div class="help-content">','<a class="less-help" href="javascript:void(0);">',c("Less"),"</a>","</div>","</div>",'<div class="alert alert-dismissable">','<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>','<span class="alert-message"></span>',"</div>",'<div class="column-headers vertically-spaced flex-column-container">','<div class="forward-column flex-column column">','<div class="column-header">','<div class="column-title">','<span class="title">',c("Unpaired forward"),"</span>",'<span class="title-info unpaired-info"></span>',"</div>",'<div class="unpaired-filter forward-unpaired-filter pull-left">','<input class="search-query" placeholder="',c("Filter this list"),'" />',"</div>","</div>","</div>",'<div class="paired-column flex-column no-flex column">','<div class="column-header">','<a class="choose-filters-link" href="javascript:void(0)">',c("Choose filters"),"</a>",'<a class="clear-filters-link" href="javascript:void(0);">',c("Clear filters"),"</a><br />",'<a class="autopair-link" href="javascript:void(0);">',c("Auto-pair"),"</a>","</div>","</div>",'<div class="reverse-column flex-column column">','<div class="column-header">','<div class="column-title">','<span class="title">',c("Unpaired reverse"),"</span>",'<span class="title-info unpaired-info"></span>',"</div>",'<div class="unpaired-filter reverse-unpaired-filter pull-left">','<input class="search-query" placeholder="',c("Filter this list"),'" />',"</div>","</div>","</div>","</div>"].join("")),middle:_.template(['<div class="unpaired-columns flex-column-container scroll-container flex-row">','<div class="forward-column flex-column column">','<ol class="column-datasets"></ol>',"</div>",'<div class="paired-column flex-column no-flex column">','<ol class="column-datasets"></ol>',"</div>",'<div class="reverse-column flex-column column">','<ol class="column-datasets"></ol>',"</div>","</div>",'<div class="flexible-partition">','<div class="flexible-partition-drag" title="',c("Drag to change"),'"></div>','<div class="column-header">','<div class="column-title paired-column-title">','<span class="title"></span>',"</div>",'<a class="unpair-all-link" href="javascript:void(0);">',c("Unpair all"),"</a>","</div>","</div>",'<div class="paired-columns flex-column-container scroll-container flex-row">','<ol class="column-datasets"></ol>',"</div>"].join("")),footer:_.template(['<div class="attributes clear">','<div class="clear">','<label class="remove-extensions-prompt pull-right">',c("Remove file extensions from pair names"),"?",'<input class="remove-extensions pull-right" type="checkbox" />',"</label>","</div>",'<div class="clear">','<input class="collection-name form-control pull-right" ','placeholder="',c("Enter a name for your new list"),'" />','<div class="collection-name-prompt pull-right">',c("Name"),":</div>","</div>","</div>",'<div class="actions clear vertically-spaced">','<div class="other-options pull-left">','<button class="cancel-create btn" tabindex="-1">',c("Cancel"),"</button>",'<div class="create-other btn-group dropup">','<button class="btn btn-default dropdown-toggle" data-toggle="dropdown">',c("Create a different kind of collection"),' <span class="caret"></span>',"</button>",'<ul class="dropdown-menu" role="menu">','<li><a href="#">',c("Create a <i>single</i> pair"),"</a></li>",'<li><a href="#">',c("Create a list of <i>unpaired</i> datasets"),"</a></li>","</ul>","</div>","</div>",'<div class="main-options pull-right">','<button class="create-collection btn btn-primary">',c("Create list"),"</button>","</div>","</div>"].join("")),helpContent:_.template(["<p>",c(["Collections of paired datasets are ordered lists of dataset pairs (often forward and reverse reads). ","These collections can be passed to tools and workflows in order to have analyses done on each member of ","the entire group. This interface allows you to create a collection, choose which datasets are paired, ","and re-order the final collection."].join("")),"</p>","<p>",c(['Unpaired datasets are shown in the <i data-target=".unpaired-columns">unpaired section</i> ',"(hover over the underlined words to highlight below). ",'Paired datasets are shown in the <i data-target=".paired-columns">paired section</i>.',"<ul>To pair datasets, you can:","<li>Click a dataset in the ",'<i data-target=".unpaired-columns .forward-column .column-datasets,','.unpaired-columns .forward-column">forward column</i> ',"to select it then click a dataset in the ",'<i data-target=".unpaired-columns .reverse-column .column-datasets,','.unpaired-columns .reverse-column">reverse column</i>.',"</li>",'<li>Click one of the "Pair these datasets" buttons in the ','<i data-target=".unpaired-columns .paired-column .column-datasets,','.unpaired-columns .paired-column">middle column</i> ',"to pair the datasets in a particular row.","</li>",'<li>Click <i data-target=".autopair-link">"Auto-pair"</i> ',"to have your datasets automatically paired based on name.","</li>","</ul>"].join("")),"</p>","<p>",c(["<ul>You can filter what is shown in the unpaired sections by:","<li>Entering partial dataset names in either the ",'<i data-target=".forward-unpaired-filter input">forward filter</i> or ','<i data-target=".reverse-unpaired-filter input">reverse filter</i>.',"</li>","<li>Choosing from a list of preset filters by clicking the ",'<i data-target=".choose-filters-link">"Choose filters" link</i>.',"</li>","<li>Entering regular expressions to match dataset names. See: ",'<a href="https://developer.mozilla.org/en-US/docs/Web/JavaScript/Guide/Regular_Expressions"',' target="_blank">MDN\'s JavaScript Regular Expression Tutorial</a>. ',"Note: forward slashes (\\) are not needed.","</li>","<li>Clearing the filters by clicking the ",'<i data-target=".clear-filters-link">"Clear filters" link</i>.',"</li>","</ul>"].join("")),"</p>","<p>",c(["To unpair individual dataset pairs, click the ",'<i data-target=".unpair-btn">unpair buttons ( <span class="fa fa-unlink"></span> )</i>. ','Click the <i data-target=".unpair-all-link">"Unpair all" link</i> to unpair all pairs.'].join("")),"</p>","<p>",c(['You can include or remove the file extensions (e.g. ".fastq") from your pair names by toggling the ','<i data-target=".remove-extensions-prompt">"Remove file extensions from pair names?"</i> control.'].join("")),"</p>","<p>",c(['Once your collection is complete, enter a <i data-target=".collection-name">name</i> and ','click <i data-target=".create-collection">"Create list"</i>. ',"(Note: you do not have to pair all unpaired datasets to finish.)"].join("")),"</p>"].join(""))};(function(){jQuery.fn.extend({hoverhighlight:function l(n,m){n=n||"body";if(!this.size()){return this}$(this).each(function(){var p=$(this),o=p.data("target");if(o){p.mouseover(function(q){$(o,n).css({background:m})}).mouseout(function(q){$(o).css({background:""})})}});return this}})}());var d=function a(n,l){l=_.defaults(l||{},{datasets:n,oncancel:function(){Galaxy.modal.hide()},oncreate:function(){Galaxy.modal.hide();Galaxy.currHistoryPanel.refreshContents()}});if(!window.Galaxy||!Galaxy.modal){throw new Error("Galaxy or Galaxy.modal not found")}var m=new k(l);Galaxy.modal.show({title:"Create a collection of paired datasets",body:m.$el,width:"80%",height:"800px",closing_events:true});m.render();window.PCC=m;return m};return{PairedCollectionCreator:k,pairedCollectionCreatorModal:d}});