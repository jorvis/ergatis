// var nodes = new Object();

function AjaxQueue( max_calls ) {
    this.calls = new Array();
    this.postProcessors = new Array();
    
    if ( ! max_calls ) {
        max_calls = 0;
    }
    
    // max number of simultaneous calls pending (0 for unlimited)
    this._max_calls = max_calls;
    this._calls_out = 0;
    
    // methods
    this.addCall = add_call;
    this.addPostProcessor = add_post_processor;
    this.start = purge_queue;
    
    // accessors
    this.getCallCount = _get_call_count;
}

function add_call( fref ) {
    this.calls.push( fref );
}

function add_post_processor( fref ) {
    this.postProcessors.push( fref );
}

function _get_call_count() {
    return this.calls.length;
}

function purge_queue() {
    if ( this.calls.length > 0 ) {
        var calls_to_make = this._max_calls - this._calls_out;
        
        for ( var i=0; i<calls_to_make; i++ ) {
            // element 0 is the send call (with arguments), element 1 is the object ref
            var fref = this.calls.shift();
            
            // this will only happen if the number in the queue is less than the max simultaneous call count
            if ( typeof fref == "undefined" ) {
                continue;
            }
            
            /* I really want it to handle the call back itself via a second passed
               parameter, but couldn't get it to work */
            /*
            // remember the original callback
            var old_callback = arr[1].onreadystatechange;
            
            // create a new one
            arr[1].onreadystatechange = function() {
                this._calls_out--;
                old_callback();
            };
            */
            
            // fire the original call
            this._calls_out++;
            fref();
        }
        
        // wait some interval, then attempt to purge again
        //  'this' is special, and we can't use it directly here or it will get
        //  evaluated as the window Obj.
        var thisObj = this;
        window.setTimeout( function(){thisObj.start()}, 1000);
        
    } else {
        // no calls left, execute any post processors
        for ( i in this.postProcessors ) {
            this.postProcessors[i]();
        }
    }
}
