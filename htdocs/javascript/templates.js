addEvent( window, 'load', initializeTemplatePage, false );

var color_invalid = 'rgb(255,150,150)';

function initializeTemplatePage() {
    addEvent( getObject('project'), 'change', selectProject, false );
    getObject('project').selectedIndex = 0;
}

function selectPipelineTemplate(path) {
    var sel = getObject('project');
    var project = sel.options[ sel.selectedIndex ].value;

    if ( checkInstructions() ) {
        getObject('autoload_template').value = path;
        document.project_selection.submit();
    }
}

function selectProject() {
    var sel = getObject('project');
    var project = sel.options[ sel.selectedIndex ].value;
    
    if ( isEmpty(project) ) {
        // hide the project templates
        getObject('project_pipelines').style.display = 'none';
        getObject('no_project').style.display = 'block';
        
    } else {
        getObject('project').style.backgroundColor = 'rgb(255,255,255)';
        
        // pull the project template list
        
        function ajaxBindCallback() {
            // progressive transitions are from 0 .. 4
            if (ajaxRequest.readyState == 4) {
                // 200 is the successful response code
                if (ajaxRequest.status == 200) {
                    ajaxCallback ( ajaxRequest.responseText );
                } else {
                    // error handling here
                    alert("there was a problem getting the template lists");
                }
            }
        }

        var ajaxRequest = null;
        var ajaxCallback = populateProjectTemplates;
        var url = './get_pipeline_templates.cgi?path=' + project + '/workflow/project_saved_templates';

        // bind the call back, then do the request
        if (window.XMLHttpRequest) {
            // mozilla, firefox, etc will get here
            ajaxRequest = new XMLHttpRequest();
            ajaxRequest.onreadystatechange = ajaxBindCallback;
            ajaxRequest.open("GET", url , true);
            ajaxRequest.send( null );

        } else if (window.ActiveXObject) {
            // IE, of course, has its own way
            ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");

            if (ajaxRequest) {
                ajaxRequest.onreadystatechange = ajaxBindCallback;
                ajaxRequest.open("GET", url , true);
                ajaxRequest.send( );
            }
        }
    }
}

function populateProjectTemplates( template_html ) {
    // populate it
    getObject('project_pipelines').innerHTML = template_html;
    
    // hide some things
    getObject('no_project').style.display = 'none';
    
    // show the templates
    if ( isEmpty(template_html) ) {
        getObject('no_project_templates_found').style.display = 'block';
    } else {
        getObject('no_project_templates_found').style.display = 'none';
        getObject('project_pipelines').style.display = 'block';
    }
}

function checkInstructions() {
    var sel = getObject('project');
    var project = sel.options[ sel.selectedIndex ].value;
    
    if ( isEmpty( project ) ) {
        // highlight the instructions
        getObject('instructions').style.fontWeight = 'bold';
        getObject('project').style.backgroundColor = color_invalid;
        getObject('instructions').focus();
        
        return false;
        
    } else {
        getObject('instructions').style.fontWeight = 'normal';
        getObject('project').style.backgroundColor = 'rgb(255,255,255)';
        return true;
    }
}
