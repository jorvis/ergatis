function showDialogue(dialogue_name) {
    hideAllDialogues();
    getObject(dialogue_name).style.display = 'block';
}

function hideAllDialogues() {
    hideDialogue('archive_verification_dialogue');
    hideDialogue('delete_dialogue');
    hideDialogue('temp_archive_verification_dialogue');
}

function hideDialogue(dialogue_name) {
    getObject(dialogue_name).style.display = 'none';
}
