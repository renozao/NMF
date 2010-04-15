
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

function acknowledgment(){
	?>
	<h3>Acknowledgements</h3>
	<div id="support">
	<b>Funding:</b>
	This work was funded by the South-African National Bioinformatics Network.
	Cathal Seoighe is funded through Science Foundation Ireland (Grant number 07/SK/M1211b).
	</div>
	<table width="100%" cellpadding="0" cellspacing="0" align="center">
		<tr>
			<td><img src="nbn.png" alt="" /></td>
			<td><a href="http://www.sfi.ie" target="_new_nmf"><img src="sfi.png" alt="" /></a><br /></td>		
		
		<td><a href="http://cbio.uct.ac.za" target="_new_nmf"><img src="cbio.png" alt="" /></a></td>
		<td><a href="http://www.uct.ac.za" target="_new_nmf"><img src="uct.gif" alt="" /></a></td>
		</tr>
	</table>
	<?
}

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

define(PKG_ROOT, '../pkg/');
define(PKG_SVN_VIEW_ROOT, 'http://r-forge.r-project.org/plugins/scmsvn/viewcvs.php/*checkout*/pkg/');

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title>
	<?php // echo $group_name; ?>
	NMF: a flexible R package for Nonnegative Matrix Factorization
	</title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php
/*
if ($handle=@fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } 
*/ 
?>

<!-- end of project description -->

<style type="text/css">

body{
	text-align: center;
	max-width:950px;
	min-width:950px;
	margin-left:auto;
	margin-right:auto;
}

a img{
	border:none;
}

div#wrapper{
	text-align: left;
}

div.h2{
	margin:19px 0px 19px 0px;
	font-size:16pt;
	font-weight:bold;
	background-color:#135CAE;
	color:white;
	padding:10px;
}

div.h3{
	margin:19px 0px 19px 0px;
	font-size:12pt;
	font-weight:bold;
}

h3, div.h3{
	background-color:#dddddd;
	padding:2px;
}

div#abstract{
	font-size: 10pt;
}

div#support{
	font-size: 10pt;
}

img{
	vertical-align: middle;
}

.comment{
	color: #770000;
}

div.menu th{
	padding:8px;
	border-top:1px solid gray;
	border-bottom:1px solid gray;
}

th.link{
	padding:8px;
	border-top:1px solid gray;
	border-bottom:1px solid gray;
	border-left:1px dotted gray;
	background-color:#dddddd;
}
div.menu th.link a{
	color:#111199;
}

div.menu td{
	border-bottom:1px solid gray;
	border-left:1px solid gray;
}

div.menu th:first-child {
	border-left:1px solid gray;
	background-color:#dddddd;
} 

div.menu th.link:hover, div.menu th.link:hover a{
	background-color:#135CAE;
	color: white;
}

div.menu th.selected, div.menu th.selected a{
	background-color:#ffffff;
	color: black;
	border-bottom:none;
}

.textbf{
	font-weight: bold;
}


</style>

<div id="wrapper">
<div class="h2">NMF: A flexible R package for Nonnegative Matrix Factorization</div>
<table width="100%" cellpadding="0" cellspacing="0">
<tr>
<td valign="top">
Renaud Gaujoux&sup1;, Cathal Seoighe*&sup2;
<br />
<small>
<p>(1) Institute of Infectious Disease and Molecular Medicine, University of Cape Town, Observatory 7925, South Africa<br />
(2) School of Mathematics, Statistics and Applied Mathematics, National University of Ireland Galway, Ireland<br />
* Corresponding author
</p>
<p>Email: Renaud Gaujoux - &lt;renaud at cbio.uct.ac.za&gt;; Cathal Seoighe - &lt;cathal.seoighe at nuigalway.ie&gt;</p>
</small>
</td>
<td><img src="nmf.png" alt="" /></td>
</tr></table>

<?php
$view = $_GET['view']; 
if( !in_array($view,array('soft', 'references')) ) 
	$view = 'soft';
?>
<div class='menu'>
<table cellpadding="0" cellspacing="0"><tr>
<th>&nbsp;</th>
<!-- <th class="link <?php echo !$view ? "selected": ""; ?>"><a href="?">Home</a></th> -->
<th class="link <?php echo $view=='soft' ? "selected": ""; ?>"><a href="?view=soft">Software</a></th>
<th class="link <?php echo $view=='references' ? "selected": ""; ?>"><a href="?view=references">References</a></th>
<td width="100%"></td>
</tr></table>
</div>
<?php

switch( $view ){
		
	case 'soft':
	?>
	<table width="100%"><tr>
<td valign="top" width="38%">
<h3>Software &amp; Documentations</h3>

<?php
function get_local_version($os='', $name = '', $version='', $rforge=true){

	if( !$os ){ // return last version for all os
		$os_array = array('nix', 'win', 'mac');
		$res = array();
		foreach( $os_array as $os)
			$res[$os] = get_local_version($os, $name, $version, $rforge);
		return $res;
	}

	// retrieve the current DESCRIPTION file
	if( !$name || !$version ){
		$desc = file_get_contents(PKG_ROOT.'DESCRIPTION');	
		// get version
		$pattern = "/Version *: *([0-9.]+)/m";
		preg_match_all($pattern, $desc, $matches);
		$version = $matches[1][0];	
		// get package name
		$pattern = "/^Package *: *([^ \n]+)/";
		preg_match_all($pattern, $desc, $matches);
		$name = $matches[1][0];	
	}
	
	$pkg_name = $name."_".$version;
	echo $rforge;			
	if( $rforge ){
		global $domain;	
		$url = 'http://'.$domain;
		if( $os == 'nix' ) $url .= "/src/contrib/".$pkg_name.".tar.gz";
		else if( $os == 'win' ) $url .= "/bin/windows/contrib/latest/".$pkg_name.".zip";
		else if( $os == 'mac' ) $url .= "/bin/macosx/leopard/contrib/latest/".$pkg_name.".tgz";
		return $url;	
	}else{		
		$url = "";
		if( $os == 'nix' ) $url = "release/".$pkg_name.".tar.gz";
		else if( $os == 'win' ) $url = "release/".$pkg_name.".zip";
		else if( $os == 'mac' ) $url = "release/".$pkg_name.".tgz";
		return $url;
	}
		
}

$local_pkgs = get_local_version('', 'NMF', '0.4', false);
function link_package($file, $name=''){
	$name = ( $name ? $name : basename($file) );
	if( !file_exists($file) ) return "<i>$name</i> [not built yet]";
	return "<a href=\"$file\">$name</a>";
}
?>
<ul>
<!-- Not working links...
<li>Package source: <a href="<?php echo $local_pkgs['nix'];?>"><?php echo basename($local_pkgs['nix']);?></a></li>
<li>MacOS X binary: <a href="<?php echo $local_pkgs['mac'];?>"><?php echo basename($local_pkgs['mac']);?></a></li>
<li>Windows binary: <a href="<?php echo $local_pkgs['win'];?>"><?php echo basename($local_pkgs['win']);?></a></li>
-->
<li>Package source: <?php echo link_package($local_pkgs['nix']);?></a></li>
<li>MacOS X binary: <?php echo link_package($local_pkgs['mac']);?></li>
<li>Windows binary: <?php echo link_package($local_pkgs['win']);?>
<br />[R devel: <a href="<?php echo "devel/i386/".basename($local_pkgs['win']);?>">i386</a> - <a href="<?php echo "devel/i386/".basename($local_pkgs['win']);?>">x86_64</a>]</li>
<li>Reference manual: <a href="NMF-manual.pdf">NMF-manual.pdf</a></li>
<li>Vignette: <a href="<?php echo PKG_SVN_VIEW_ROOT?>inst/doc/NMF-vignette.pdf?root=nmf">NMF-vignette.pdf</a></li>
<li>News/ChangeLog:	<a href="<?php echo PKG_SVN_VIEW_ROOT?>NEWS?root=nmf">NEWS</a></li>
</ul>
<table><tr>
<td><a href="http://www.r-project.org" target="_new_nmf"><img src="Rlogo.jpg" alt="" /></a></td>
<td><a href="http://www.bioconductor.org" target="_new_nmf"><img src="bioconductor.png" alt="" /></a></td>
</tr>
</table>
</td>
<td valign="top" width="62%">
<?php
ini_set('allow_url_fopen', 'on');
function http_get($query)
{
	// in case we run on localhost (not CBIO)
	//if( $_SERVER['HTTP_HOST'] == 'localhost' ) return file_get_contents($query);
	
	require_once "HTTP/Request.php";
	
	$req =& new HTTP_Request($query);
	//$req->setProxy("campusnet.uct.ac.za", 8080, "gjxren001", "kngprw10#");
	
	// check for error
	if ( PEAR::isError($req->sendRequest()) ) 
	{    	
		echo "Could not retrieve query: $query";	
		return;
	}
	
	$res = $req->getResponseBody();
	return $res;
}

function get_cran_version($pkg){
	
	$url = "http://cran.r-project.org/package=$pkg";
	$old = ini_set('default_socket_timeout', 5);
	$cran = http_get($url);
	ini_set('default_socket_timeout', $old);
	if( $cran && ereg("$pkg_([0-9.]+)\.tar.gz", $cran, $version) )
		return $version[1];
}
//$latest_version = get_cran_version('NMF');
$latest_version = '0.3.3';
?>
<h3>Install&nbsp;&amp;&nbsp;Updates
<font style="font-size:10pt;color:#990000"><?php echo $latest_version ? "[CRAN version: $latest_version]" : "" ?></font>
</h3>
<p>The NMF package is still under active development. The latest version of the package is available on CRAN, please click <a target="_new_nmf" href="http://cran.r-project.org/package=NMF">here</a> for a direct link to the dedicated web page on the central mirror.</p>
<p>Installation is done via the following standard call in the R console:</p>
<code><pre style="background-color:#aaccaa;padding:10px;border:1px solid gray">
<font class="comment"># install the NMF package from default CRAN mirror</font>
install.packages('NMF');
</pre></code>

<h3 style="margin-bottom:0px">Project summary</h3>
<!-- R-Forge Logo -->
<p>
The NMF package has a dedicated <strong>project summary page</strong> on <a target="_new" href="http://r-forge.r-project.org">R-Forge</a>.
You can find it <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>.
</p>
<p>It provides access to useful tools such as the source code SVN repository, trackers to enter bugs or request features, and a mailing-list.</p>

</td>
</tr></table>

	<?php
	
		// Add acknowledgment
		acknowledgment();
		break;
		
	//END_SOFTWARE
	case 'references':
		readfile('references.html');
		break;
	default:
	if( false ):
	?>
	
<h3>Abstract</h3>
<div id="abstract">
<p>
<b>Background:</b>
Nonnegative Matrix Factorization (NMF) is an unsupervised learning technique that has been applied successfully in several fields, including signal processing, facial recognition and text mining. 
Recent applications of NMF in bioinformatics have demonstrated its ability to extract meaningful information from high-dimensional data such as gene expression microarrays. Developments in NMF theory and applications have resulted in a variety of algorithms and methods. 
However, most NMF implementations have been on commercial platforms, while those that are freely available typically require programming skills.
This limits their use by the wider research community.
</p>

<p>
<b>Results:</b>
Our objective is to provide the bioinformatics community with an open-source, easy-to-use and unified interface to standard NMF algorithms, as well as with a simple framework to help implement and test new NMF methods.
For that purpose, we have developed a package for the R/BioConductor platform. The package ports public code to R, and is structured to enable users to easily modify and/or add algorithms.
It includes a number of published NMF algorithms and initialization methods and facilitates the combination of these to produce new NMF strategies.
Commonly used benchmark data and visualization methods are provided to help in the comparison and interpretation of the results.
</p>

<p>
<b>Conclusions:</b>
The NMF package helps realize the potential of Nonnegative Matrix Factorization, especially in bioinformatics, providing easy access to methods that have already yielded new insights in many applications.
Documentation, source code and sample data are available from CRAN at <a target="_new_nmf" href="http://cran.r-project.org">http://cran.r-project.org</a>.
</p>
</div>

<?php 
// Add acknowledgment
acknowledgment();
endif;

} //end switch
?>
</div>

</body>
</html>
