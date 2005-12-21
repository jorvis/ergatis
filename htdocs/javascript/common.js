/*
    this document is for common script portions and shouldn't contain any
    completely ergatis-specific methods.  most should be common house-keeping
    and utility functions.
*/

// are we dealing with a browser that needs hand-holding?
var IE = (navigator.appName == "Microsoft Internet Explorer") ? true : false;


// returns an object reference by its id
function getObject(obj) {
    // the W3C way
    if (document.getElementById) {
        return document.getElementById(obj);
    
    // the netscape way
    } else if (document.layers) {
        return eval("document." + obj);
    }
}
