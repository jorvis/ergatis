$(document).ready(function() {
    $(document).keypress(function(e) {
        if (e.keyCode == 13) {
            document.login_form.submit();
        }
    });     
});
