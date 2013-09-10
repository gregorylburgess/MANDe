
Ncss.map = null;

Ncss.log("Ncss loading...");

Ncss.initMapPreview = function(){

	var vector_layer = new OpenLayers.Layer.Vector({
		renderers:["Canvas"],
	});	


	var wkt = new OpenLayers.Format.WKT();
	vector_layer.addFeatures( wkt.read(gridWKT) );

	var geometry  = OpenLayers.Geometry.fromWKT(gridWKT);
	var gridCentroid  = geometry.getCentroid();

	//Checks if we got a geometry, if so we should get a good centroid for it.
	//otherwise preview is not available
	//if(typeof gridCentroid !== "undefined" ){
	if(gridCentroid !== null ){

//		Init map preview
		map = new OpenLayers.Map({
			div : "gridPreview",
			theme : null,
			controls :[new OpenLayers.Control.Navigation(), new OpenLayers.Control.Zoom(), new OpenLayers.Control.MousePosition({numDigits:2, separator: '|', emptyString:'Mouse is not over map', formatOutput:ncssFormatOutput} )], //No controls in preview
			layers : [
			          new OpenLayers.Layer.MapServer("Basic", "http://vmap0.tiles.osgeo.org/wms/vmap0", 
			        		  {layers:'basic'},
			        		  {wrapDateLine : true}
			          )
			          ],		
			          center: new OpenLayers.LonLat(gridCentroid.x,gridCentroid.y),
			          zoom:0
		});	

		map.addLayer(vector_layer);

		Ncss.gridCrossesDateLine = ( vector_layer.getDataExtent().left <= 180 & vector_layer.getDataExtent().right >= 180);

		console.log("layer data extent:"+ vector_layer.getDataExtent().left+", "+vector_layer.getDataExtent().bottom+", "+vector_layer.getDataExtent().right+", "+vector_layer.getDataExtent().top );

		map.zoomToExtent( new OpenLayers.Bounds( vector_layer.getDataExtent().left, vector_layer.getDataExtent().bottom, vector_layer.getDataExtent().right, vector_layer.getDataExtent().top ) );

		
	}else{
		//preview not available --> hide the map div
		$('#gridPreviewFrame').css('display','none');
	}

};


 Ncss.rezoom = function(){	
	alert("moving to: " + $('[name="west"]').val() + " " + $('[name="east"]').val());
	map.addLayer(vector_layer);
	map.zoomToExtent( new OpenLayers.Bounds($('[name="west"]').val(),
											$('[name="south"]').val(),
											$('[name="east"]').val(),
											$('[name="north"]').val() 
											)
					);		
}

var ncssFormatOutput = function(lonLat) {
    var digits = parseInt(this.numDigits);
    
    if(Ncss.gridCrossesDateLine){
    	if(lonLat.lon < 0 ){
    		lonLat.lon+=360;
    	}
    }	
    var newHtml =
        this.prefix +
        lonLat.lon.toFixed(digits)+
        this.separator + 
        lonLat.lat.toFixed(digits) +
        this.suffix;
    return newHtml;
}
