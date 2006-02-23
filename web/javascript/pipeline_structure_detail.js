var notes_visible = 0;

function hideShow(ident) {
    if ( document.getElementById(ident).style.display == 'block' ){
        document.getElementById(ident).style.display = 'none';
        document.images["i" + ident].src = "../images/plus.png";
    } else {
        document.getElementById(ident).style.display = 'block';
        document.images["i" + ident].src = "../images/minus.png";
    }
}

function toggleNotes() {

    divs = document.getElementsByTagName("div");

    for (i=0; i < divs.length; i++) {
        // skip this div if it isn't the right class
        if (divs[i].getAttribute("class") == 'changenote') {
            if (notes_visible) {
                divs[i].style.display = 'none';
            } else {
                divs[i].style.display = 'block';
            }
        }

    }

    if (notes_visible) {
        notes_visible = 0;
    } else {
        notes_visible = 1;
    }
}
