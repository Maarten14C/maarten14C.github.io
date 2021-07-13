function popupURL(newurl)
{
// var newurl = popupURL.arguments;

   var availableWidth = screen.availWidth - 100;
   var availableHeight = screen.availHeight - 100;

   windowLeft = 50;
   windowTop = 50;

   newWindow = window.open(newurl,"newwin"
                 ,"width=" + availableWidth 
              + ",height=" + availableHeight 
              + ",left=" + windowLeft 
              + ",top=" + windowTop 
              + ",location=yes,resizable=yes,scrollbars=yes,toolbar=yes");
    newWindow.focus();

}

function MM_goToURL() 
{ //v3.0
  var i, args=MM_goToURL.arguments;
  document.MM_returnValue = false;

  for(i=0; i<(args.length-1); i+=2)
     eval(args[i]+".location='"+args[i+1]+"'");
}

function MM_findObj(n, d) 
{ //v4.01
   var p,i,x;
   if(!d) d=document;
   if((p=n.indexOf("?"))>0&&parent.frames.length)
   {
       d=parent.frames[n.substring(p+1)].document;
       n=n.substring(0,p);
   }
   if(!(x=d[n])&&d.all)
      x=d.all[n];
   for(i=0;!x&&i<d.forms.length;i++)
      x=d.forms[i][n];

   for(i=0;!x&&d.layers&&i<d.layers.length;i++)
      x=MM_findObj(n,d.layers[i].document);

   if(!x && d.getElementById)
      x=d.getElementById(n);
   return x;
}

function MM_showHideLayers() 
{ //v6.0

  var i,p,v,obj,args=MM_showHideLayers.arguments;
  for (i=0; i<(args.length-2); i+=3)
     if((obj=MM_findObj(args[i]))!=null)
     {
        v=args[i+2];
        if(obj.style) 
        {
           obj=obj.style; 
           v=(v=='show')?'visible':(v=='hide')?'hidden':v;
        }
        obj.visibility=v;
     }
}
