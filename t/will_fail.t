use Math::Pari ':all';

print "1..4\n";

$x = PARIvar 'x';
$y = PARIvar 'y';

# In fact this would segfault if moved to the end...
$bnf = buchinitfu($x**2-$x-57,0.2,0.2);
print "# bnf=$bnf\n";
print "# bnf[6]=$bnf->[6]\n";

eval {
  isideal($bnf->[-1+7],PARImat_tr([[5,1], [0,1]]));
};
$test++;
chomp($@), print "# err=$@\nnot " if $@;
print "ok $test\n";

$qpol=$y**3-$y-1;
setrand(1);
$bnf2=buchinit($qpol);
$nf2=$bnf2 -> [-1+7];

$un=mod(1,$qpol);
$w=mod($y,$qpol);
$p=$un*($x**5-5*$x+$w);

$aa = rnfpseudobasis($nf2,$p);
$test++;
$aaa = pari_print($aa);
$bbb = '[[[1, 0, 0]~, [0, 0, 0]~, [0, 0, 0]~, [-2, 0, 0]~, [11, 0, 0]~; [0, 0, 0]~, [1, 0, 0]~, [0, 0, 0]~, [2, 0, 0]~, [-8, 0, 0]~; [0, 0, 0]~, [0, 0, 0]~, [1, 0, 0]~, [1, 0, 0]~, [4, 0, 0]~; [0, 0, 0]~, [0, 0, 0]~, [0, 0, 0]~, [1, 0, 0]~, [-2, 0, 0]~; [0, 0, 0]~, [0, 0, 0]~, [0, 0, 0]~, [0, 0, 0]~, [1, 0, 0]~], [[1, 0, 0; 0, 1, 0; 0, 0, 1], [1, 0, 0; 0, 1, 0; 0, 0, 1], [1, 0, 0; 0, 1, 0; 0, 0, 1], [1, 0, 3/5; 0, 1, 2/5; 0, 0, 1/5], [1, 0, 8/25; 0, 1, 22/25; 0, 0, 1/25]], [416134375, 212940625, 388649575; 0, 3125, 550; 0, 0, 25], [-1280, 5, 5]~]';
print "# aaa=$aaa\n# exp=$bbb\nnot " unless $aaa eq $bbb;
print "ok $test\n";

# If uncommented, it segfaults
$test++;
# print "# Skipping next one, would segfault:\nnot ok $test\n";
eval { $d = rnfdiscf($nf2,$p) };
print "not " if $@ or $d ne PARI '[[416134375,212940625,388649575;0,3125,550;0,0,25],[-1280,5,5]~]';
print "ok $test\n";

$nfpol = $x**5 - 5 * $x**3 + 5 * $x + 25;
$nf = initalg($nfpol);
$vp = primedec($nf,3)->[-1+1];
$ans = idealhermite2($nf,$vp->[-1+2],3);
$aaa = pari_print $ans;
$bbb = '[3, 1, 2, 2, 2; 0, 1, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 1, 0; 0, 0, 0, 0, 1]';

$test++;
print "# aaa=$aaa\n# exp=$bbb\nnot " unless $aaa eq $bbb;
print "ok $test\n";

print STDERR "Expect that most tests here fail.  Fix them if you can...\n";
