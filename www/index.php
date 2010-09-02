<?php

session_start();

// force check on CRAN
//$_SESSION['cran_cache'] = 0;
if( isset($_GET['check_cran']) ){	
	$_SESSION['cran_cache'] = 0;
	header('location: ?view=soft');
}

define('LAST_VERSION', '0.4.6');
define('LAST_CRAN_VERSION', '0.4.6');
define('CRAN_MIRROR', "http://cran.r-project.org");
define('PKGNAME', "NMF");

?>
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

define('PKG_ROOT', '../pkg/');
define('PKG_SVN_VIEW_ROOT', 'https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/');

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
	<link href="<?php echo $themeroot; ?>css/theme.css" rel="stylesheet" type="text/css" />
	<script type="text/javascript" src="jquery.js"></script>  
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
	margin:19px 0px 10px 0px;
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

div#citation{
	font-size:12pt;	
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
<div id="citation"><i>BMC Bioinformatics</i> 2010, <b>11</b>:367 -
<a href="http://www.biomedcentral.com/1471-2105/11/367" target="_new_nmf">doi:10.1186/1471-2105-11-367</a> - PMCID: PMC2912887
<br />
<br />
</div>
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
<p>Authors: Renaud Gaujoux - &lt;renaud at cbio dot uct dot ac dot za&gt;; Cathal Seoighe - &lt;cathal dot seoighe at nuigalway dot ie&gt;</p>
</small>
</td>
<td><img src="nmf.png" alt="" /></td>
</tr></table>

<?php
$view = isset($_GET['view']) ? $_GET['view'] : ''; 
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

ini_set('allow_url_fopen', 'on');
function http_get($query)
{
	// in case we run on localhost (not CBIO)
	//if( $_SERVER['HTTP_HOST'] == 'localhost' ) return file_get_contents($query);
	return file_get_contents($query);
	
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

function get_cran_version($pkg, $default){
	
	$url = CRAN_MIRROR."/package=$pkg";
	$old = ini_set('default_socket_timeout', 5);
	$cran = http_get($url);
	//$cran_xml = simplexml_load_string($cran);
	//echo $cran_xml->body;		
	preg_match_all("/<p>(.*)<\/p>/siU", $cran, $desc);
	preg_match_all("/<table .*<\/table>/siU", $cran, $downloads);
	$pkg_info = array('desc' => $desc[0][0],'detail' => $downloads[0][0], 'links' => $downloads[0][1]);	
	ini_set('default_socket_timeout', $old);
	$pkg_version = ( $cran && ereg($pkg."_([0-9.]+)\.tar\.gz", $cran, $version) ? $version[1] : $default);
	return array($pkg_version, $pkg_info);
}
	
if( !(!isset($_SESSION['CRAN']) || $cran_info=$_SESSION['CRAN']) || (!isset($_SESSION['cran_cache']) || $_SESSION['cran_cache'] + (60*5) < time()) ){
	$cran_info = get_cran_version('NMF', LAST_CRAN_VERSION);	
	$_SESSION['CRAN'] = $cran_info;
	$_SESSION['cran_cache'] = time();
}
$latest_version = $cran_info[0];

//print_r($cran_info);
//print the package description
echo $cran_info[1]['desc'];
?>

<table width="100%"><tr>
<td valign="top" width="38%">

<h3>Software &amp; Documentations</h3>

<?php

// change table into ul
$download_links = preg_replace("`<tr[^>]*> *<td[^>]*>(.*)</td> *<td[^>]*>(.*)</td> *</tr>`siU", '<li>\1 \2</li>', $cran_info[1]['links']);
// fix links
$download_links = preg_replace('`href=[\'"](\.[^\'"]+)[\'"]`siU','href="'.CRAN_MIRROR.'/web/packages/'.PKGNAME.'/\1"',$download_links);
echo "<ul>$download_links</ul>";
?>

<a href="#" onclick="$('#pkg_details').slideToggle();">+ Details</a><br /><br />
<div id="pkg_details" style="display:none"><?php 
$details_links = preg_replace('`href=[\'"](\.[^\'"]+)[\'"]`siU','href="'.CRAN_MIRROR.'/web/packages/'.PKGNAME.'/\1"',$cran_info[1]['detail']);
echo $details_links;?></div>
<br />

<?php
/* OLD CODE FOR LINKS

function get_local_version($os='', $name = '', $default_version='', $rforge=true){

	if( !$os ){ // return last version for all os
		$os_array = array('nix', 'win', 'mac');
		$res = array();
		foreach( $os_array as $os)
			$res[$os] = get_local_version($os, $name, $default_version, $rforge);
		return $res;
	}

	// find the last local version
	$tars = scandir('release');
	$maxv = '';
	foreach( $tars as $t ){
		if( !eregi("^[a-z.]+_([0-9.\-]+)\.tar\.gz$", $t, $v) ) continue;
		if( version_compare($v[1], $maxv) > 0 ) 
			$maxv = $v[1];
	}
	if( $maxv ) 
		$version = $maxv;
	else 
		$version = $default_version;
		
	// retrieve the current DESCRIPTION file
	if( !$name || !$version ){
		//$desc = file_get_contents(PKG_ROOT.'DESCRIPTION');
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

$local_pkgs = get_local_version('', 'NMF', LAST_VERSION, false);
function link_package($file, $name='', $alt='[not built yet]'){
	$name = ( $name ? $name : basename($file) );
	$alt = ' '.$alt;
	if( !file_exists($file) ) return "<i>$name</i>$alt";
	return "<a href=\"$file\">$name</a>";
}

 * ?>
<ul>
<!-- Not working links...
<li>Package source: <a href="<?php echo $local_pkgs['nix'];?>"><?php echo basename($local_pkgs['nix']);?></a></li>
<li>MacOS X binary: <a href="<?php echo $local_pkgs['mac'];?>"><?php echo basename($local_pkgs['mac']);?></a></li>
<li>Windows binary: <a href="<?php echo $local_pkgs['win'];?>"><?php echo basename($local_pkgs['win']);?></a></li>
-->
<li>Package source: <?php echo link_package($local_pkgs['nix']);?></a></li>
<li>MacOS X binary: <?php echo link_package($local_pkgs['mac']);?></li>
<li>Windows binary: <?php echo link_package("release/i386/".basename($local_pkgs['win']), "i386",'[N/A]');?> - <?php echo link_package("release/x86_64/".basename($local_pkgs['win']), "x86_64", '[N/A]');?>
<br />[R devel: <?php echo link_package("devel/i386/".basename($local_pkgs['win']), "i386",'[N/A]');?> - <?php echo link_package("devel/x86_64/".basename($local_pkgs['win']), "x86_64", '[N/A]');?></li>
<li>Reference manual: <a href="NMF-manual.pdf">NMF-manual.pdf</a></li>
<li>Vignette: <a href="<?php echo PKG_SVN_VIEW_ROOT?>inst/doc/NMF-vignette.pdf?root=nmf">NMF-vignette.pdf</a></li>
<li>News/ChangeLog:	<a href="<?php echo PKG_SVN_VIEW_ROOT?>NEWS?root=nmf">NEWS</a></li>
</ul>
*/
?>
<table><tr>
<td><a href="http://www.r-project.org" target="_new_nmf"><img src="Rlogo.jpg" alt="" /></a></td>
<td><a href="http://www.bioconductor.org" target="_new_nmf"><img src="bioconductor.png" alt="" /></a></td>
</tr>
</table>
</td>

<td valign="top" width="62%">
<h3>Install&nbsp;&amp;&nbsp;Updates
<font style="font-size:10pt;color:#990000">
<?php echo $latest_version ? "[CRAN version: <span id=\"version\">$latest_version</span>]" : "" ?>
&nbsp;<a style="font-size:8pt" href="?check_cran=1" onclick="document.getElementById('version').innerHTML='.....';">Recheck</a>
</font>
</h3>
<p>The NMF package is still under active development. The latest stable version of the package is available on CRAN.
 <a target="_new_nmf" href="http://cran.r-project.org/package=NMF">Go to the <b>CRAN</b> page</a>.</p>
<p>Installation is done via the following standard call in the R console:</p>
<code><pre style="background-color:#aaccaa;padding:10px;border:1px solid gray">
<font class="comment"># install the NMF package from default CRAN mirror</font>
install.packages('NMF');
</pre></code>

<h3 style="margin-bottom:0px">Project & Development</h3>
<!-- R-Forge Logo -->
<p>
The NMF package has a dedicated <strong>project on <a target="_new" href="http://r-forge.r-project.org">R-Forge</a></strong>.
 <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/">Go to the <strong>project home</strong></a>.
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
