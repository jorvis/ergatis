/*
   Written by Terry Friesen,  tfriesen@mts.net
   http://www.mts.net/~tfriesen/dhtml/
   
   this script gives Netscape 6 (and mozilla browsers) the following IE methods:
        removeNode()
        replaceNode()
        swapNode()
        applyElement()
        contains()
        insertAdjacentText()
        insertAdjacentHTML()
            examples: 
            foo.insertAdjacentHTML('BeforeBegin', '<bar>');
            foo.insertAdjacentHTML('AfterBegin', '<bar>');
            foo.insertAdjacentHTML('BeforeEnd', '<bar>');
        insertAdjacentElement()
*/

if(self.Node&&self.Node.prototype){
    Node.prototype.removeNode=remove_Node;
    Node.prototype.replaceNode=replace_Node;
    Node.prototype.swapNode=swap_Node;
    Element.prototype.applyElement=apply_Element;
    Element.prototype.contains=_contains;
    Element.prototype.insertAdjacentText=insertAdj_Text;
    Element.prototype.insertAdjacentHTML=insertAdj_HTML;
    Element.prototype.insertAdjacentElement=insertAdj_El;
    Element.prototype.insert__Adj=insert__Adj;
}

function remove_Node(a1){
var p=this.parentNode;
if(p&&!a1){
var df=document.createDocumentFragment();
for(var a=0;a<this.childNodes.length;a++){
df.appendChild(this.childNodes[a])
}
p.insertBefore(df,this)
}
return p?p.removeChild(this):this;
}

function replace_Node(a1) {
    return this.parentNode.replaceChild(a1,this);
}

function swap_Node(a1){
var p=a1.parentNode;
var s=a1.nextSibling;
this.parentNode.replaceChild(a1,this);
p.insertBefore(this,s)
return this;
}

function apply_Element(a1,a2){
if(!a1.splitText){
a1.removeNode();
if(a2&&a2.toLowerCase()=="inside"){
for(var a=0;a<this.childNodes.length;a++){
a1.appendChild(this.childNodes[a])
}
this.appendChild(a1)
}
else{
var p=this.parentNode;
p.insertBefore(a1,this);
a1.appendChild(this);
}
return a1;
}
}

function _contains(a1){
var r=document.createRange();
r.selectNode(this);
return r.compareNode(a1)==3;
}

function insertAdj_Text(a1,a2){
var t=document.createTextNode(a2||"")
this.insert__Adj(a1,t);
}

function insertAdj_HTML(a1,a2){
var r=document.createRange();
r.selectNode(this);
var t=r.createContextualFragment(a2);
this.insert__Adj(a1,t);
}

function insertAdj_El(a1,a2){
this.insert__Adj(a1,a2);
return a2;
}

function insert__Adj(a1,a2){
var p=this.parentNode;
var s=a1.toLowerCase();
if(s=="beforebegin"){p.insertBefore(a2,this)}
if(s=="afterend"){p.insertBefore(a2,this.nextSibling)}
if(s=="afterbegin"){this.insertBefore(a2,this.childNodes[0])}
if(s=="beforeend"){this.appendChild(a2)}
}
