#! perl
use Math::Pari ':all';

eval { link_gnuplot() };
if ($@ =~ m%^Can't locate Term/Gnuplot.pm in \@INC%) {
  print STDERR "# Can't locate Term/Gnuplot.pm in \@INC, ignoring the test\n";
  print "1..0\n";
  exit;
} elsif ($@) {
  die $@;
} else {
  print "1..1\n";
}
setprecision 9;
$x = PARIvar 'x';

$t = plothsizes();
die if Math::Pari::typ($t) < 17;
$w=floor($t->[0]*0.42)-1;
$h=floor($t->[1]*0.42)-1;
$dw=floor($t->[0]*0.05)+1;
$dh=floor($t->[1]*0.05)+1;
plotinit(2, 2*$w+10, 2*$h+10);
plotinit(3, $w, $h);
plotrecth(3, $x, -5, 5, sub {sin($x)}, 2, 0);
plotcopy(3, 2, $dw, $dh);
plotinit(3, $w, $h);
plotrecth(3, $x, -5, 5, sub {[sin($x),cos(2*$x)]}, 0, 0);
plotcopy(3, 2, $w + 2*$dw, $dh);
plotinit(3, $w, $h);
plotrecth(3, $x, -5, 5, sub {[sin(3*$x), cos(2*$x)]}, 1, 0);
plotcopy(3, 2, $dw, $h + 2*$dh);
plotinit(3, $w, $h);
plotrecth(3, $x, -5, 5, sub {[sin($x), cos($x), sin(3*$x),cos(2*$x)]}, 1, 0);
plotcopy(3, 2, $w+2*$dw, $h+2*$dh);
plotdraw([2, 0, 0]);

print "ok 1\n";
print STDERR "Press ENTER\n";
<>;
