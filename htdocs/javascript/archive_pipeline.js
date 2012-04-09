var settimmer = 0;

$(document).ready(function() {
    // Our pipeline list will be stored in a hidden input
    var pipelineListURL = $("#pipeline_list_input").val();

    var wintimer = window.setInterval(function() {
	    var timeCounter = $("b[id=redirect_count]").html();
	    var updateTime = eval(timeCounter)- eval(1);
	    $("b[id=redirect_count]").html(updateTime);

	    if(updateTime <= 0){
            window.clearInterval(wintimer);
            window.location = pipelineListURL;
	    }
    }, 1000);
});
