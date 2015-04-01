
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<! --- R-Forge Logo --- >
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<h1>Packages</h1>

<h2>marelac: Tools for Aquatic Sciences</h2>

<p>Datasets, constants, conversion factors, utilities for MArine,
Riverine, Estuarine, LAcustrine and Coastal science.</p>

<p>The package contains among others:</p>

<ol>

<li> chemical and physical constants and datasets, e.g. atomic
weights, gas constants, the earths bathymetry;

<li> conversion factors (e.g. gram to mol to liter, barometric units,
temperature, salinity),</li>

<li> physical functions, e.g. to estimate concentrations of
conservative substances, gas transfer and diffusion coefficients, the
Coriolis force and gravity,</li>

<li>thermophysical properties of the seawater, as from the UNESCO
polynomial or from the more recent derivation based on a Gibbs
function.</li>

</ol>

<h2>marelacTeaching</h2>

<ul>
<li> includes lecture notes "Using R for scientific computing" for the novice non-statistician R-user. 
</ul>

<h2>Links</h2>

<ul>
<li> <a href="http://cran.r-project.org/web/packages/marelac/"><strong>Package description on CRAN</strong></a>. </p>
<li> <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>Project summary page on R-Forge</strong></a>. </p>
</ul>

</body>
</html>
