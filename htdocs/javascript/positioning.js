/*
    this document contains methods for detecting the position of object such
    as the mouse or document elements.
*/

// follow the mouse
if (! IE) {
    document.captureEvents(Event.MOUSEMOVE);
}

// set the tracking of the mouse
var mouseX = 0;
var mouseY = 0;
document.onmousemove = storeMouseXY;

// gets the current horizontal position of the mouse
function getMouseX(e) {
    if (IE) {
        return event.clientX + document.body.scrollLeft;
    } else {
        return e.pageX;
    }
}


// gets the current vertical position of the mouse
function getMouseY(e) {
    if (IE) {
        return event.clientY + document.body.scrollTop;
    } else {
        return e.pageY;
    }
}

// returns the horizontal position of a passed element, in pixels
function getObjectX(e) {
    var posX = e.offsetLeft;
    
    // assumes you're not asking for the position of the whole document (duh)
    e = e.offsetParent;
    
    // we have to build the position by traversing up the DOM.
    while (e != null) {
        posX += e.offsetLeft;
        e = e.offsetParent;
    }
    
    return posX;
}

// returns the vertical position of a passed element, in pixels
function getObjectY(e) {
    var posY = e.offsetTop;
    
    // assumes you're not asking for the position of the whole document (duh)
    e = e.offsetParent;
    
    // we have to build the position by traversing up the DOM.
    while (e != null) {
        posY += e.offsetTop;
        e = e.offsetParent;
    }
    
    return posY;
}

// gets the width of the browser window - 20 (a sad attempt to adjust for the
//   scroll bar.)   Returns 300 if it fails to get the size (because surely the
//   browser is at least that wide.
function getWindowWidth() {
    // there are three ways to do this
    
    // mozilla way
    if ( document.body.offsetLeft ) {
        return document.body.offsetLeft - 20;
        
    // IE way
    } else if ( document.body.offsetWidth ) {
        return document.body.offsetWidth - 20;
    
    // netscape way
    } else if ( window.innerWidth ) {
        return window.innerWidth - 20;
    
    } else {
        // what fallback should we use here?
        alert("error, unable to determine screen center.  don't panic.");
        return 300;
    }
}

// gets the height of the browser window - 20 (a sad attempt to adjust for the
//   scroll bar.)   Returns 300 if it fails to get the size (because surely the
//   browser is at least that wide.
function getWindowHeight() {
    // there are three ways to do this
    
    // mozilla way
    if ( document.body.offsetTop ) {
        return document.body.offsetTop - 20;
    
    // IE way
    } else if ( document.body.offsetHeight ) {
        return document.body.offsetHeight - 20;
    
    // netscape way
    } else if ( window.innerHeight ) {
        return window.innerHeight - 20;
    
    } else {
        // what fallback should we use here?
        alert("error, unable to determine screen center.  don't panic.");
        return 300;
    }
}

// sets the X and Y coordinates of the mouse to the mouseX and mouseY globals
// returns true for the fun of it
function storeMouseXY(e) {
    mouseX = getMouseX(e);
    mouseY = getMouseY(e);
    
    // for debugging
    // document.title = "(" + mouseX + ", " + mouseY + ")";
    return true;
}

// sets the posX attribute of the object.
function storeObjectX(e) {
    e.posX = getObjectX(e);
}

// sets the posY attribute of the object.
function storeObjectY(e) {
    e.posY = getObjectY(e);
}

// sets the posX and posY attributes of the object.
function storeObjectXY(e) {
    storeObjectX(e);
    storeObjectY(e);
}




















