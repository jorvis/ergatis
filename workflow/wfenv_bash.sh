WF_WORK=/usr/local/devel/ANNOTATION/workflow
export WF_WORK

WF_TEMPLATE=$WF_ROOT/templates
export WF_TEMPLATE

WF_JARS=$WF_ROOT/jars
export WF_JARS

cp=$WF_JARS/wf.jar:$WF_JARS/castor-0.9.4.1-xml.jar:$WF_JARS/castor-0.9.4.1.jar:$WF_JARS/xerces.jar:$WF_JARS/jconn2.jar
export cp

PATH=${WF_ROOT}:${WF_ROOT}/bin:${WF_ROOT}/add-ons/bin:/usr/local/java/1.4.1/bin/:$PATH
export PATH

