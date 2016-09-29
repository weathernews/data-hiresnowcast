#! /usr/bin/env perl
use File::Basename;
use Geo::Proj4;
use GPX;

#
# Google Maps Tile projection parameters
#

$proj = Geo::Proj4->new(proj   => "merc",
			ellps  => "sphere",
			lat_0  => 0,
			lon_0  => 0);
# tile size
($imgw, $imgh) = (256,256);
# world 
$g_n_lat	= 85.05;
$g_w_lon	= -180;
$g_s_lat	= -85.05;
$g_e_lon	= 180;

my($x1, $y1) = $proj->forward($g_n_lat, $g_w_lon);
my($x2, $y2) = $proj->forward($g_s_lat, $g_w_lon);
my($x3, $y3) = $proj->forward($g_s_lat, $g_e_lon);
my($x4, $y4) = $proj->forward($g_n_lat, $g_e_lon);

$g_lx = ($x1 < $x2) ? $x1 : $x2;
$g_rx = ($x3 > $x4) ? $x3 : $x4;
$g_ty = ($y1 > $y4) ? $y1 : $y4;
$g_by = ($y2 < $y3) ? $y2 : $y3;


$orgw = 159;
$orgh = 159;
$o_n_lat = 25.998958;	# one of some block
$o_w_lon = 123.001563;
$o_s_lat = 25.667708;
$o_e_lon = 123.498438;

#$orgw = 2560;
#$orgh = 3360;
$jo_n_lat =  48;	# JMA radar
$jo_w_lon = 118;
$jo_s_lat =  20;
$jo_e_lon = 150;


# （緯度方向は等間隔ではないので）１タイルあたりの経度幅を求める
$o_hreso = ($o_e_lon - $o_w_lon) / $orgw * $imgw;
# そのタイルより細かくなる解像度 (Level of Detail) を求める
$max_lod = int(log(360 / $o_hreso) / log(2)) + 1;
$ntile = 2 ** $max_lod;


@llpts = ();
open(LL,"10-grib-latlon.txt") || die $!;
while (<LL>) {
    my ($lat0,$lon0,$lat1,$lon1) = split(/[, \r\n]+/,$_);
    $lat0 /= 1000000;
    $lon0 /= 1000000;
    $lat1 /= 1000000;
    $lon1 /= 1000000;
    push(@llpts, [$lat0,$lon0]);
    push(@llpts, [$lat1,$lon1]);
    push(@llpts, [$lat0,$lon1]);
    push(@llpts, [$lat1,$lon0]);
}
close(LL);

open(O,">tile_info.h");
print O << "+++";
typedef struct _tile {
  char id[16];
  int lt_lat,lt_lon,rb_lat,rb_lon;
  unsigned char *p;
} TILE;
TILE tile[] = {
+++
    ;

$n = 0;
for ($i = 0; $i < $ntile; $i++) {
    for ($j = 0; $j < $ntile; $j++) {
	($lt_lat, $lt_lon, $rb_lat, $rb_lon) = tile_boundingbox($max_lod,$i,$j);

	if (($rb_lat > $jo_n_lat) || 
	    ($lt_lat < $jo_s_lat) ||
	    ($rb_lon < $jo_w_lon) ||
	    ($lt_lon > $jo_e_lon)) {
	    next;
	}

	$found = 0;
	foreach $llpt (@llpts) {
	    my ($lat,$lon) = ($llpt->[0], $llpt->[1]);
	    if (($lt_lat > $lat) && ($lat > $rb_lat) &&
		($lt_lon < $lon) && ($lon < $rb_lon)) {
		$found = 1;
		last;
	    }
	}
	next if ($found == 0);

	$fn = sprintf("%d_%d_%d",$max_lod,$i,$j);

	$lat0 = int($lt_lat * 1000000);
	$lon0 = int($lt_lon * 1000000);
	$lat1 = int($rb_lat * 1000000);
	$lon1 = int($rb_lon * 1000000);
	
	print O qq(${delim}{"$fn",$lat0,$lon0,$lat1,$lon1});
	$delim = ",\n";
	$n++;
    }
}
print O qq(};\n);
print O qq(const int tile_count = $n;\n);
close(O);
exit;


sub tile_boundingbox {
    my($lv,$tilex,$tiley) = @_;
    my $dv = 2 ** $lv;
    if (($tilex >= $dv) || ($tiley >= $dv)) {
	print STDERR "tile_boundingbox : addressing error\n";
	return;
    }
    my $tilew = ($g_rx - $g_lx) / $dv;
    my $tileh = ($g_ty - $g_by) / $dv;

    my $x0 = $tilew * $tilex + $g_lx;
    my $y0 = $g_ty - $tileh * $tiley;
    my $x1 = $x0 + $tilew;
    my $y1 = $y0 - $tileh;

    my ($lt_lat, $lt_lon) = $proj->inverse($x0, $y0);
    my ($rb_lat, $rb_lon) = $proj->inverse($x1, $y1);
    return ($lt_lat, $lt_lon, $rb_lat, $rb_lon);
}

