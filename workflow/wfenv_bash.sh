WF_WORK=/usr/local/projects/PANDA/workflow
export WF_WORK

WF_TEMPLATE=$WF_ROOT/templates
export WF_TEMPLATE

WF_JARS=$WF_ROOT/jars
export WF_JARS

cp=$WF_JARS/wf.jar:$WF_JARS/castor-0.9.4.1-xml.jar:$WF_JARS/castor-0.9.4.1.jar:$WF_JARS/xerces.jar:$WF_JARS/jconn2.jar
export cp

