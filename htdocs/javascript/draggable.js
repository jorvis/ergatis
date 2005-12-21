// required scripts:
//      positioning.js

var drag = false;
var elemX = 0;
var elemX = 0;

function cancelDrag(e) {
    e.style.zIndex -= 10;

    if (drag) {
        window.clearInterval(drag);
    }
}

function startDrag(e) {
    // bring the object above other things on the page.
    e.style.zIndex += 10;
    theID = e.id;

    var fix = storeObjectXY(e);
    difX = e.posX - mouseX;
    difY = e.posY - mouseY;
    
    if (drag) {
        window.clearInterval(drag);
    }

    if (! IE) {
        drag = setInterval("dragObject(getObject(theID))");
    } else {
        drag = setInterval("dragObject(" + theID + ")", 1);
    }
}

function dragObject(e) {
    e.style.left = mouseX + difX;
    e.style.top = mouseY + difY;

    // drag any children

    // check to see if it is within the drop zone of any element
}
