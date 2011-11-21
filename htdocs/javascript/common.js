/*
    allows us to easily center any element selected using jquery
*/
jQuery.fn.center = function () {
    this.css("position","absolute");
    this.css("top", (($(window).height() - this.outerHeight()) / 2) + $(window).scrollTop() + "px");
    this.css("left", (($(window).width() - this.outerWidth()) / 2) + $(window).scrollLeft() + "px");
    return this;
}

/*
    code allowing spin.js to interface with jquery
*/
$.fn.spin = function(opts) {
  this.each(function() {
    var $this = $(this),
        data = $this.data();

    if (data.spinner) {
      data.spinner.stop();
      delete data.spinner;
    }
    if (opts !== false) {
      data.spinner = new Spinner($.extend({color: $this.css('color')}, opts)).spin(this);
    }
  });
  return this;
};

/*
* Makes an input element look un-editable
*/
function inputDisplayUneditable(element) {
    jQuery(element).css('background-color', 'rgb(255, 255, 255)');               
    jQuery(element).css('borderTop', '0px');               
    jQuery(element).css('borderLeft', '0px');               
    jQuery(element).css('borderRight', '0px');               
    jQuery(element).css('borderBottom', '0px');               

}

/*
* Makes an input element look editable
*/
function inputDisplayEditable(element) {
    jQuery(element).css('background-color', 'rgb(230,230,230)');               
    jQuery(element).css('borderTop', '1px solid rgb(160,160,160)');               
    jQuery(element).css('borderLeft', '1px solid rgb(160,160,160)');               
    jQuery(element).css('borderRight', '1px solid rgb(160,160,160)');               
    jQuery(element).css('borderBottom', '1px solid rgb(160,160,160)');               
}

/*
    this document is for common script portions and shouldn't contain any
    completely ergatis-specific methods.  most should be common house-keeping
    and utility functions or functions javascript should have that don't.
    
    a lot of the code is mine, some isn't.  i've tried to include information
    about other others whenever available in the comments above each function.
*/

// are we dealing with an unruly browser that needs hand-holding?
var IE = (navigator.appName == "Microsoft Internet Explorer") ? true : false;

/*
    cross-browser event handling for IE5+, NS6+ and Mozilla/Gecko
    by Scott Andrew
    http://www.scottandrew.com/weblog/articles/cbs-events
*/
function addEvent( elm, evType, fn, useCapture ) {
    if ( elm.addEventListener ) {
        elm.addEventListener( evType, fn, useCapture );
        return true;
        
    } else if ( elm.attachEvent ) {
        var r = elm.attachEvent( 'on' + evType, fn );
        return r;
        
    } else {
        elm[ 'on' + evType ] = fn;
    }
}

/*
    clearInput simply clears the value of an input element inputId if, optionally, its 
    current value is clearVal.  this  allows you to have text boxes with initial values like
    'enter text here', which you can pass to this function to clear upon first entry.
    
    sets the element's value to empty string.
    
    - jorvis
*/
function clearInput( inputId, clearVal ) {
    if ( isEmpty(clearVal) ) {
        getObject(inputId).value = '';
        
    } else if ( getObject(inputId).value == clearVal ) {
        getObject(inputId).value = '';
    }
}

/*
    getElementsByClassName written by Jonathan Snook, http://www.snook.ca/jonathan
    Add-ons by Robert Nyman, http://www.robertnyman.com
    
    Some ways to call it

    To get all a elements in the document with a "info-links" class.
        getElementsByClassName(document, "a", "info-links");
    
    To get all div elements within the element named "container", with a "col" class.
        getElementsByClassName(document.getElementById("container"), "div", "col"); 
    
    To get all elements within in the document with a "click-me" class.
        getElementsByClassName(document, "*", "click-me"); 
*/

function getElementsByClassName(oElm, strTagName, strClassName){
    var arrElements = (strTagName == "*" && oElm.all)? oElm.all : oElm.getElementsByTagName(strTagName);
    var arrReturnElements = new Array();
    strClassName = strClassName.replace(/\-/g, "\\-");
    var oRegExp = new RegExp("(^|\\s)" + strClassName + "(\\s|$)");
    var oElement;
    for(var i=0; i<arrElements.length; i++){
        oElement = arrElements[i];      
        if(oRegExp.test(oElement.className)){
            arrReturnElements.push(oElement);
        }   
    }
    return (arrReturnElements)
}

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

// checks the many ways that a field can be empty
function isEmpty( val ) {
    var rgx = new RegExp('^\\s+$');

    // check for null
    // check for empty string
    // check for all whitespace
    if ( val == null ||
         val == undefined ||  
         val == '' || 
         val.match(rgx) ) {
        
        return true;    
    } else {
        return false;
    }
}

function debug( msg ) {
    getObject('debug_box').innerHTML += msg + '<br>';
}


/*
 * Copyright 2005 Matthew Eernisse (mde@fleegix.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Original code by Matthew Eernisse (mde@fleegix.org), March 2005
 * Additional bugfixes by Mark Pruett (mark.pruett@comcast.net), 12th July 2005
 * Multi-select added by Craig Anderson (craig@sitepoint.com), 24th August 2006
 *
 * Version 1.3
*/

/**
 * Serializes the data from all the inputs in a Web form
 * into a query-string style string.
 * @param docForm -- Reference to a DOM node of the form element
 * @param formatOpts -- JS object of options for how to format
 * the return string. Supported options:
 *    collapseMulti: (Boolean) take values from elements that
 *    can return multiple values (multi-select, checkbox groups)
 *    and collapse into a single, comman-delimited value
 *    (e.g., thisVar=asdf,qwer,zxcv)
 * @returns query-string style String of variable-value pairs
 */
function formData2QueryString(docForm, formatOpts) {
  
  var opts = formatOpts || {};
  var str = '';
  var formElem;
  var lastElemName = '';
  
  for (i = 0; i < docForm.elements.length; i++) {
    formElem = docForm.elements[i];
    
    switch (formElem.type) {
      // Text fields, hidden form elements
      case 'text':
      case 'hidden':
      case 'password':
      case 'textarea':
      case 'select-one':
        str += formElem.name + '=' + encodeURIComponent(formElem.value) + '&'
        break;
        
      // Multi-option select
      case 'select-multiple':
        var isSet = false;
        for(var j = 0; j < formElem.options.length; j++) {
          var currOpt = formElem.options[j];
          if(currOpt.selected) {
            if (opts.collapseMulti) {
              if (isSet) {
                str += ',' + encodeURIComponent(currOpt.value);
              }
              else {
                str += formElem.name + '=' + encodeURIComponent(currOpt.value);
                isSet = true;
              }
            }
            else {
              str += formElem.name + '=' + encodeURIComponent(currOpt.value) + '&';
            }
          }
        }
        if (opts.collapseMulti) {
          str += '&';
        }
        break;
      
      // Radio buttons
      case 'radio':
        if (formElem.checked) {
          str += formElem.name + '=' + encodeURIComponent(formElem.value) + '&'
        }
        break;
        
      // Checkboxes
      case 'checkbox':
        if (formElem.checked) {
          // Collapse multi-select into comma-separated list
          if (opts.collapseMulti && (formElem.name == lastElemName)) {
            // Strip of end ampersand if there is one
            if (str.lastIndexOf('&') == str.length-1) {
              str = str.substr(0, str.length - 1);
            }
            // Append value as comma-delimited string
            str += ',' + encodeURIComponent(formElem.value);
          }
          else {
            str += formElem.name + '=' + encodeURIComponent(formElem.value);
          }
          str += '&';
          lastElemName = formElem.name;
        }
        break;
        
    }
  }
  // Remove trailing separator
  str = str.substr(0, str.length - 1);
  return str;
}
