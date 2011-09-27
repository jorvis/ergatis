/*
    assumes the page has a 
        textarea with id pipeline_comment
*/
var pcomment;

function initializeCommentContainer () {
    pcomment = getObject('comment_shadow').innerHTML;

    if ( isEmpty( pcomment ) ) {
        getObject('comment_add').style.display = 'inline';
    
    } else {
        getObject('comment_add').style.display = 'none';
        getObject('comment_shadow').style.display = 'block';
    }
    
    addEvent( getObject('comment_add'), 'click', editComment, false );
    addEvent( getObject('comment_cancel'), 'click', cancelComment, false );
    addEvent( getObject('comment_save'), 'click', saveComment, false );
    addEvent( getObject('comment_shadow'), 'click', editComment, false );
    addEvent( getObject('comment_shadow'), 'mouseover', lookEditable, false );
    addEvent( getObject('comment_shadow'), 'mouseout', lookUneditable, false );
}

function cancelComment() {
    // hide edit box
    getObject('comment').style.display = 'none';

    // display the shadow
    getObject('comment_shadow').style.display = 'block';
    
    if ( isEmpty( pcomment ) ) {
        getObject('comment_add').style.display = 'inline';
    
    } else {
        getObject('comment_add').style.display = 'none';
        getObject('comment_shadow').style.display = 'block';
    }
    
    // hide the save and cancel labels
    getObject('comment_save').style.display = 'none';
    getObject('comment_cancel').style.display = 'none';
    
}

function editComment() {
    // populate the edit box
    getObject('comment').value = pcomment;

    // hide the shadow
    getObject('comment_shadow').style.display = 'none';

    // display edit box
    getObject('comment').style.display = 'block';
    
    // toss cursor
    getObject('comment').focus();
    
    // hide click label
    getObject('comment_add').style.display = 'none';
    
    // display save/cancel labels
    getObject('comment_save').style.display = 'inline';
    getObject('comment_cancel').style.display = 'inline';
}

function lookEditable() {
    var container = getObject('comment_shadow');
    
    container.style.backgroundColor = 'rgb(230,230,230)';
    container.style.borderTop = '1px solid rgb(160,160,160)';
    container.style.borderLeft = '1px solid rgb(160,160,160)';
    container.style.borderRight = '1px solid rgb(200,200,200)';
    container.style.borderBottom = '1px solid rgb(200,200,200)';
}

function lookUneditable() {
    var container = getObject('comment_shadow');
    
    container.style.backgroundColor = 'rgb(255,255,255)';
    container.style.border = '1px solid rgb(255,255,255)';
}

function saveComment() {
    getObject('user_message_container').innerHTML = 'saving pipeline comment';
    getObject('user_message_container').style.display = 'block';

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
                ajaxCallback ( ajaxRequest.responseText );
            } else {
                // error handling here
                alert("there was a problem saving the pipeline comment");
            }
        }
    }

    var ajaxRequest = null;
    var ajaxCallback = _save_comment_handler;
    var url = './save_comment.cgi';
    var form_string = formData2QueryString( document.forms['comment_form'] );

    // bind the call back, then do the request
    if (window.XMLHttpRequest) {
        // mozilla, firefox, etc will get here
        ajaxRequest = new XMLHttpRequest();
        ajaxRequest.onreadystatechange = ajaxBindCallback;
        ajaxRequest.open("POST", url , true);
        ajaxRequest.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
        ajaxRequest.setRequestHeader("Connection", "close");
        ajaxRequest.send( form_string );
        
    } else if (window.ActiveXObject) {
        // IE, of course, has its own way
        ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");

        if (ajaxRequest) {
            ajaxRequest.onreadystatechange = ajaxBindCallback;
            ajaxRequest.open("POST", url , true);
            ajaxRequest.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
            ajaxRequest.setRequestHeader("Connection", "close");
            ajaxRequest.send( form_string );
        }
    }

    // display 'pipeline build progress saved' label
    getObject('user_message_container').innerHTML = 'comment saved';
    window.setTimeout( "getObject('user_message_container').style.display = 'none'", 2000);
}

function _save_comment_handler() {
    pcomment = getObject('comment').value;
    getObject('comment_shadow').innerHTML = pcomment;
    
    cancelComment();
}

/*
    fun with listeners
*/
addEvent( window, 'load', initializeCommentContainer, false );













