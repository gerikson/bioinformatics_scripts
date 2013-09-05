<%@ page language="java" contentType="text/html; charset=ISO-8859-1"
    pageEncoding="ISO-8859-1"%>
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">

<html>

<head>
<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">
<title>Scripps Genome ADVISER | Annotation and Distributed Variant Interpretation Server</title>
<LINK REL=STYLESHEET HREF="all_styles.css" TYPE="text/css">
</head>
<body>
<div id="content-wrapper">
<%@ include file="top_navigation.jsp" %>

<script src="http://www.genomics.scripps.edu/ADVISER/scripts/gatag.js" type="text/javascript"></script>
<script type="text/javascript">

  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-39253100-1']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();

</script>

<div id="page-content">
<br/>
<h1>Privacy Tool</h1>
<p style="margin:18px;">
The Privacy Tool extracts the genotype(s) from the VCF/CGI files and 
implants clinically associated variants into the list of transmitted variants.
For the clinically associated variants we used UCSC table named Flagged SNPs132.  
The Privacy Tool accepts VCF files.  
</br></br>
Download the <a href="Apps/dif_privacy.zip"  target="_blank" onclick="_gaq.push(['_trackEvent','Privacy Tool download','ZIP','Privacy Tool download']);">Privacy Tool</a>.
</br></br>
Questions, Comments, Suggestions? Post on <a href="Forum.jsp"> Biostar</a>.
</p>
<p style="margin:18px;">
A deidentified file ready for SG ADVISER annotation will be *Your 
Filename*.Deidentified.tar.gz. To annotate this file you can proceed to genomics.scripps.edu  
and upload the file for annotation.
</p>
<p style="margin:18px;">
Identification strips the clinically associated variants previously implanted and imports 
ack genotype(s). For identification please import the SG ADVISER annotated file then from 
the tool folder select  the *filename*.Identified file. The final annotated and identified
file will be  final_*filename*.Identified. You can import this file into the SG ADVISER UI
for sorting and filtering.
</p>
<p style="margin:18px;">
NOTE: depending of the file(s) size both process can take a long time (up to 10 minutes). 
Please be patient and wait for a notification that deidentification/ identification has 
finished and you see the output files in the tool folder.
</p> 
</body>
</html>
</div>
<%@ include file="footer.jsp" %>
</div>
</body>

</html>
