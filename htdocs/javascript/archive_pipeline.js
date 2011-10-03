window.onload = function() {
    // Our pipeline list will be stored in a hidden input
    var pipelineListURL = document.getElementById('pipeline_list_input').value;

    setTimeout(function() {
        window.location.href=pipelineListURL;
    }, 5000);
}
