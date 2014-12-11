define(["mvc/list/list-item","mvc/dataset/dataset-list","mvc/base-mvc","utils/localization"],function(d,f,b,c){var e=d.FoldoutListItemView;var a=e.extend({className:e.prototype.className+" job",id:function(){return["job",this.model.get("id")].join("-")},foldoutPanelClass:f.DatasetList,initialize:function(g){if(g.logger){this.logger=this.model.logger=g.logger}this.log(this+".initialize:",g);e.prototype.initialize.call(this,g);this.tool=g.tool||{};this.jobData=g.jobData||{};this.linkTarget=g.linkTarget||"_blank"},_swapNewRender:function(g){e.prototype._swapNewRender.call(this,g);if(this.model.has("state")){this.$el.addClass("state-"+this.model.get("state"))}return this.$el},_getFoldoutPanelOptions:function(){var g=e.prototype._getFoldoutPanelOptions.call(this);return _.extend(g,{collection:this.model.outputCollection,selecting:false})},_labelParamMap:function(){var h=this.model.get("params"),g={};_.each(this.tool.inputs,function(j){if(j.label&&j.model_class!=="DataToolParameter"){g[j.label]=h[j.name]}});return g},_labelInputMap:function(){var g=this,h={};_.each(this.jobData.inputs,function(i){var j=g._findToolInput(i.name);if(j){h[j.label]=i}});return h},_findToolInput:function(g){var i=this.tool.inputs,h=_.findWhere(i,{name:g});if(h){return h}return this._findRepeatToolInput(g,i)},_findRepeatToolInput:function(i,j){j=j||this.tool.inputs;var g=_.find(j,function(k){return i.indexOf(k.name)===0});if(!g){return undefined}var h=_.find(g.inputs,function(k){return i.indexOf(k.name)!==-1});return h},toString:function(){return"JobListItemView("+this.model+")"}});a.prototype.templates=(function(){var g=b.wrapTemplate(['<div class="list-element">','<div class="id"><%= model.id %></div>','<div class="warnings"></div>','<div class="selector">','<span class="fa fa-2x fa-square-o"></span>',"</div>",'<div class="primary-actions"></div>','<div class="title-bar"></div>','<div class="details"></div>',"</div>"]);var i=b.wrapTemplate(['<div class="title-bar clear" tabindex="0">','<div class="title">','<span class="name"><%- view.tool.name %></span>',"</div>",'<div class="subtitle">','<span class="description"><%- view.tool.description %></span','<span class="create-time">'," ",c("Created"),": <%= new Date( job.create_time ).toString() %>, ","</span","</div>","</div>"],"job");var j=b.wrapTemplate(['<div class="subtitle">','<span class="description"><%- view.tool.description %></span',"</div>"],"job");var h=b.wrapTemplate(['<div class="details">','<div class="params">',"<% _.each( view._labelInputMap(), function( input, label ){ %>",'<div class="input" data-input-name="<%= input.name %>" data-input-id="<%= input.id %>">','<label class="prompt"><%= label %></label>','<span class="value"><%= input.content.name %></span>',"</div>","<% }) %>","<% _.each( view._labelParamMap(), function( param, label ){ %>",'<div class="param" data-input-name="<%= param.name %>">','<label class="prompt"><%= label %></label>','<span class="value"><%= param %></span>',"</div>","<% }) %>","</div>","</div>"],"job");return _.extend({},e.prototype.templates,{titleBar:i,subtitle:j,details:h})}());return{JobListItemView:a}});