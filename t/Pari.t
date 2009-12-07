#!/usr/bin/perl

BEGIN { unshift @INC, './lib', '../lib';
    require Config; import Config;
}

$test = 0;
$| = 1;
print "1..",&last,"\n";

sub test {
  $test++; if (shift) {print "ok $test\n";1} else {print "not ok $test\n";0}
}

use Math::Pari;

test(1);			# 1

$a=new Math::Pari 1;

test($a);			# 2
test(ref $a eq "Math::Pari");	# 3

{
  package Math::Pari;
  ::test(defined &pari2nv); # 4
}

test(defined &Math::Pari::pari2nv); # 5
test(Math::Pari::pari2iv($a)==1); # 6

# This looks like a bug in the Pari conversion to double:

test(abs(Math::Pari::pari2nv($a)-1.0)<2**-49); # 7
#####print "`", Math::Pari::pari2nv($a)-1.0, "'\n";

$add = Math::Pari::loadPari("gadd");

test($add);			# 8
test("$add" =~ /^CODE/);	# 9

$b= &$add (34,67);

test($b);			# 10
test(ref $b eq "Math::Pari");	# 11
test(Math::Pari::pari2iv($b)==101); # 12

$b=Math::Pari::gadd(11,12);	# Uninitialized value here

test($b);			# 13
test(ref $b eq "Math::Pari");	# 14
test(Math::Pari::pari2iv($b)==23); # 15

$c = Math::Pari::gsub(17,93);

test($c);			# 16
test(ref $c eq "Math::Pari");	# 17
test(Math::Pari::pari2iv($c)==-76); # 18

$c = do { package Math::Pari; gsub(19,93) }; # And here (No MORE :-)

test($c);			# 19
test(ref $c eq "Math::Pari");	# 20
test(Math::Pari::pari2iv($c)==-74); # 21

$d=Math::Pari::gadd($b,$c);

test(ref $d eq "Math::Pari");	# 22
test(Math::Pari::pari2iv($d)==-51); # 23

$e=new Math::Pari "matrix(2,2,j,k,j+k)";

test(ref $e eq "Math::Pari");	# 24

$cc=Math::Pari::pari2iv($c);

test(! ref $cc);		# 25

$ee=Math::Pari::pari2pv($e);

test(! ref $ee) || print ref $ee;		# 26

test($ee eq "[2,3;3,4]"); # 27

#$f=Math::Pari::matinvr($e);
$f= $e ** -1;

test(ref $f eq "Math::Pari");	# 28
test(Math::Pari::pari2pv($f) eq "[-4,3;3,-2]"); # 29

# Overload

test("$f" eq "[-4,3;3,-2]") || print "$f\n";	# 30

$f+=1;

test("$f" eq "[-3,3;3,-1]") || print "$f\n";	# 31
test($f == "[-3,3;3,-1]");	# 32

$g=(new Math::Pari "[1,2;3,2]");
$gg=Math::Pari::gpui($g,-1);

test($gg);			# 33
$g=(new Math::Pari "[1,2;3,2]")**-1;

test($g == $gg);		# 34

# Should work eventually

test("$g" eq "[-1/2,1/2;3/4,-1/4]") || print "$g\n"; # 35
#test(1);

test($g == "[-1/2,1/2;3/4,-1/4]"); # 36
#test(1);

$f  = new Math::Pari "[2,3;3,1]";
$ff = PARImat [[2,3],[3,1]];

test(Math::Pari::geq($f,$ff));	# 37
test($f==$ff);			# 38
test(Math::Pari::gne($f,$ff)==0);	# 39
test(!($f!=$ff));		# 40
test("$f" eq "$ff");		# 41
test($f eq $ff);		# 42

$f  = new Math::Pari "[2,3;3,2]";

test($f != $ff);		# 43
test(!($f == $ff));		# 44
test("$f" ne "$ff");		# 45
test($f ne $ff);		# 46

$f = new Math::Pari "[2,3;2,1]";
$fff=$f**2;

test($fff eq "[10,9;6,7]");	# 47
test($fff eq $f*$f);		# 48

test(PARI(3)<4);		# 49
test(3<PARI(4));		# 50
test(!(PARI(3)<2));		# 51
test(!(4<PARI(3)));		# 52

test(PARI(3) lt 4);		# 53
test(3 lt PARI(4));		# 54
test(!(PARI(3) lt 2));		# 55
test(!(4 lt PARI(3)));		# 56

eval <<'EOE';
use Math::Pari qw(I Pi);

#sub Math::Pari::i;
#sub Math::Pari::pi;

# Bug in Pari: depends on the order of execution

# Stupid test, probably version- or machine-specific

test(abs(exp((Pi()+7e-29)*I)+1) >= 1e-29); # 57
test(abs(exp(Pi()*I)+1) <= 1e-28); # 58

test(('x'*I)**2 == "-x^2");	# 59
EOE
;

if ( ($tot, $now, $onstack, $offstack) = Math::Pari::memUsage) {
  print "# total = $tot, now = $now, onstack = $onstack, offstack = $offstack\n";
}

sub incr1 {
  my $y = shift;
  print "# Adding $y to $x\n";
  $x += $y;
  #warn "Got $y, Sum $x\n";
}

Math::Pari::installPerlFunction(\&incr1, "incrJ");

$x = 0;
Math::Pari::fordiv(28, 'j', 'incrJ(j)');

test($x == 56) or print "$x\n";	# 60

$x = 0;
Math::Pari::fordiv(28, 'j', 'incr1(j)'); # autoloading Perl->Pari

test($x == 56);			# 61

$x = 0;
$j = PARI 'j';
Math::Pari::fordiv(28, $j, 'incr1(j)'); # autoloading Perl->Pari

test($x == 56);			# 62

$x = 0;
Math::Pari::fordiv(28, 'j', sub { incr1(PARI 'j') } ); # Perl code ->Pari expr

test($x == 56);			# 63
print "# x = '$x'\n" if $x != 56;

$x = 0;
$j = PARIvar 'j';
Math::Pari::fordiv(28, 'j', sub { incr1($j) } ); # Perl code ->Pari expr

test($x == 56);			# 64
print "# x = '$x'\n" if $x != 56;

$x = 0;
$j = PARI 'j';
Math::Pari::fordiv(28, $j, sub { incr1($j) } ); # Perl code ->Pari expr

test($x == 56);			# 65
print "# x = '$x'\n" if $x != 56;

$x = 0;
$j = 3.4;
Math::Pari::fordiv(28, $j, sub { incr1($j) } ); # Perl code ->Pari expr

test($x == 56);			# 66
print "# x = '$x'\n" if $x != 56;

$x = 0;
$j = PARI 'x+O(x^7)';
Math::Pari::fordiv(28, $j, sub { incr1($j) } ); # Perl code ->Pari expr

test($x == 56);			# 67
print "# x = '$x'\n" if $x != 56;

#use Math::Pari 'pi';

$x = 0;
$j = PARI 3.14;
Math::Pari::fordiv(28, $j, sub { incr1($j) } ); # Perl code ->Pari expr

test($x == 56);			# 68
print "# x = '$x'\n" if $x != 56;

undef $x;
eval { $x  = Pi()/0; } ;	# Needs G_KEEPERR patch!
#chomp($@), print "# '$@'\n" if $@;
print "# x = '$x'\n" if defined $x;

test( $@ =~ /^PARI:\s+\*+\s+division\s+by\s+zero\b/i ) or print "$@\n"; # 69

*O=\&Math::Pari::O;

$x = PARI 'x';
$y = 1+O($x,6);
test("$y" eq '1+O(x^6)');	# 70

$x = PARI 'x';
$y = 1+O($x**6);
test("$y" eq '1+O(x^6)');	# 71

$y = 1+O(5,6);
test("$y" eq '1+O(5^6)') or print "$y\n";	# 72

sub counter { $i += shift; }
$i = 145;
PARI 'k=5' ;
Math::Pari::fordiv(28, 'j', 'k=k+counter(j)');
test(PARI('k') == 984);		# 73

sub get2 ($$) {7}
test(PARI('get2(3,-19)==7'));	# 74

$res = 9;
eval { $res = PARI('get2(3)') };
#print "# res='$res' err=`$@'\n";
test($res == 9 and $@ =~ /expected\scharacter:\s\',/); # 75

$res = 9;
eval { $res = PARI('get2(3,1,67)') };
#print "# res='$res' err=`$@'\n";
test($res == 9 and $@ =~ /expected\scharacter:\s\'\)/); # 76

$res = Math::Pari::sum($l,5,9,sub{'kk'**$l},7);
#print "res = `$res'\n";
test("$res" eq 'kk^9+kk^8+kk^7+kk^6+kk^5+7') or print "$res\n";	# 77

$var = PARI 'var';
$var1 = PARIvar 'var';

Math::Pari::changevalue('var', 13);
test(PARI('var') eq '13');	# 78

PARI('var=12');
test(PARI('var') eq '12');	# 79

Math::Pari::changevalue($var1, 15);
test(PARI('var') eq '15');	# 80
test($var1 eq '15');		# 81

Math::Pari::changevalue($var, 51);
test(PARI('var') eq '51');	# 82

my $v = PARI [2,7,11];
$v->[1] = 17;
test("$v" eq "[2,17,11]");	# 83

my $mat = PARImat [[1,2,3],[4,8,5]];
$mat->[1][1] = 18;
test("$mat" eq "[1,4;2,18;3,5]"); # 84

$mat->[0] = [56,67,78];
test("$mat" eq "[56,4;67,18;78,5]"); # 85

eval {$mat->[0] = [56,67]};
test($@ =~ /Assignment of a columns into a matrix of incompatible height/); # 86

eval {$mat->[0] = 56};
test($@ =~ /Not a vector where column of a matrix expected/); # 87

$mat = PARImat [[1,2,3]];
eval {$mat->[0] = [56,67]};
test("$mat" eq "[56;67]");	# 88


sub last {88}
