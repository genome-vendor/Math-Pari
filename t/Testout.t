#! perl
# 	$rcs = ' $Id: testout.t,v 1.2 1997/09/22 10:13:37 ilya Exp ilya $ ' ;	

use Math::Pari qw(:DEFAULT pari_print :all);
use vars qw($x $y $z $k $t $q $a $u $j $l $name $other $n);
print "1..565\n";

push @ARGV, 'libPARI/testouts' unless @ARGV;
$| = 1;
@seen{qw(pi i euler a x y z k t q u j l n name other)} = (' ', ' ', ' ', ('$') x 100);
$x = PARI('x');
$y = PARI('y');
$z = PARI('z');
$k = PARI('k');
$t = PARI('t');
$q = PARI('q');
$a = PARI('a');
$u = PARI('u');
$j = PARI('j');
$l = PARI('l');
$n = PARI('n');
$name = PARI('name');
$other = PARI('other');

while (<>) {
  last if /^stacksize/;
}
main_loop:
while (<>) {
  s/^(\?\s+)+// or die "Malformed question: `$_'";
  next if /^\\\\/;		# Comment
  $bad = /^\\/;			# \precision = 
  $wasbadprint = /\b(plot)\b/;
  $wasprint = /\b((|p|tex)print|plot)\b/;
  chomp;
  $in = $_;
  $_ = '';
  $_ = <> while defined and /^$/; # Before warnings
  $_ = <> while defined and /^\s*\*+\s*warning/i; # skip warnings
  if (s/^\s*\*+\s*(.*)//) {		# error
    process_error($in,$_,$1);
    next;
  }
  defined or die "Can't find an answer";
  chomp;
  process_test($in, 1, []), redo if /^\?\s/; # Was a void 
  s/^%\d+\s*=\s*// or die "Malformed answer: $_" unless $bad or $wasprint;
  if ($_ eq '' or $wasprint) {	# Answer is multiline
    @ans = $_ eq '' ? () : ($_) ;
    while (<>) {
      last if /^\?\s+/;
      next if /^$/;
      chomp;
      push @ans, $_;
    }
    if ($wasbadprint) {
      process_print($in, @ans);
    } elsif ($wasprint) {
      process_test($in, 'print', [@ans]);
    } else {
      process_test($in, 0, [@ans]);
    }
    redo main_loop;
  }
  if ($bad) {
    process_set($in, $_);
  } else {
    process_test($in, 0, [$_]);
  }
}

sub format_matrix {
  my $in = shift;
  my @in = split /;/, $in;
  'PARImat_tr([[' . join('], [', @in) . ']])';
}

sub format_vvector {
  my $in = shift;
  $in =~ s/~\s*$//;
  "PARIcol($in)";
}

sub mformat {
  return join("\t", @_) unless @_ > 1 and $_[0] =~ /^\[/;
  @_ = grep {!/^$/} @_;
  return join("\t", @_) if grep {!/^\s*\[.*\]\s*$/} @_;	# Not matrix
  #return join("\t", @_) if grep {!/^\s*\([^,]*,\s*$/} @_; # Extra commas
  map {s/^\s*\[(.*)\]\s*$/$1/} @_;
  my @arr = map { join ', ', split } @_;
  '[' . join('; ', @arr) . ']';
}

sub mformat_transp {
  return join("\t", @_) unless @_ > 1 and $_[0] =~ /^\[/;
  @_ = grep {!/^$/} @_;
  return join("\t", @_) if grep {!/^\s*\[.*\]\s*$/} @_;	# Not matrix
  #return join("\t", @_) if grep {!/^\s*\([^,]*,\s*$/} @_; # Extra commas
  map {s/^\s*\[(.*)\]\s*$/$1/} @_;
  my @arr = map { [split] } @_;
  my @out;
  my @dummy = ('') x @{$arr[0]};
  for my $ind (0..$#{$arr[0]}) {
    for my $subarr (@arr) {
      @$subarr > $ind or $subarr->[$ind] = ''; 
    }
    push @out, join ', ', map {$_->[$ind]} @arr;
  }
  '[' . join('; ', @out) . ']';
}

sub massage_floats {
  my $in = shift;
  $pre = shift || "16g";
  $in =~ s/(.\d*)\s+e/$1E/gi;	# 1.74 E-78
  $in =~ s/\b(\d+\.\d*(e[-+]?\d+)?|\d{10,})\b/sprintf "%.${pre}", $1/gei;
  $in;
}

sub o_format {
  my ($var,$power) = @_;
  return " PARI('O($var^$power)') " if defined $power;
  return " PARI('O($var)') ";
}

sub process_test {
  my ($in, $noans, $out) = @_;
  my $doprint;
  $doprint = 1 if $noans eq 'print';
  $c++;
  # First a trivial processing:
  $in =~ s/\b(\d+|[a-z]+\(\))\s*\\\s*(\d+(\^\d+)?)/ gdivent($1,$2)/g; # \
  $in =~ s/\b(\d+)\s*\\\/\s*(\d+)/ gdivround($1,$2)/g; # \/
  $in =~ s/\b(\w+)\s*!/ ifact($1)/g; # !
  if ($in =~ /\\/) {		# \\ for division unsupported
    $c--;
    process_error($in, $out, '\\');
  } elsif ($in =~ /^(\w+)\s*\([^()]*\)\s*=/ and 0) { # XXXX Not implemented yet
    $c--;
    process_definition($1, $in);    
  } elsif ($in =~ /[!_\']/) {	# Factorial
    print "# Skipping (ifact/conj/deriv) `$in'\nok $c\n";
  } else {
    # work with "^", need to treat differently inside o()
    $in =~ s/\^/^^^/g;
    $in =~ s/\bo\(([^()^]*)(\^\^\^([^()]*))?\)/ o_format($1,$3) /gei;
    $in =~ s/\^\^\^/**/g;	# Now treat it outside of O()
    $in =~ s/\[([^\[\];]*;[^\[\]]*)\]/format_matrix($1)/ge; # Matrix
    $in =~ s/\[([^\[\];]*)\]\s*~/format_vvector($1)/ge; # Vertical vector
    if ($in =~ /\[[^\]]*;/) {	# Matrix
      print "# Skipping (matrix notation) `$in'\nok $c\n";
      return;
    } elsif ($in =~ /(^|[\(=])%/) {
      print "# Skipping (history notation) `$in'\nok $c\n";
      return;
    } elsif ($in =~ /
		      (
			\b 
			( while | if | for | goto | label | changevar | gettime | forstep )
			\b 
		      |
			\w+ \( \w+ \) = 
		      |
			\b 
			(
			  my _
			)?
			p? print \( 
			( \[ | (1 ,)? PARImat )
		      |
			\b forprime .* \){4}
		      )
		    /x) {
      # It is not clear why changevar gives a different answer in GP
      print "# Skipping (converting test for '$1' needs additional work) `$in'\nok $c\n";
      return;
    } elsif ($in =~ / \b ( rnfdiscf | rnfpseudobasis | idealhermite2 | isideal ) \b /x) {
#    } elsif ($in =~ / \b ( rnfhermitebasis | rnfsteinitz | idealhermite2 | isideal ) \b /x) {
      # Will result in segfault or wrong answer
      print "# Skipping (would fail, checked in different place) `$in'\nok $c\n";
      return;
    } elsif ($in =~ /\bget(heap|stack)\b/) { # Meaningless
      print "# Skipping meaningless `$in'\nok $c\n";
      return;
    } elsif ($in =~ /\b(nonesuch now)\b/) {
      print "# Skipping (possibly FATAL $1) `$in'\nok $c\n";
      return;
    }
    # Convert transposition
    $in =~ s/(\w+(\([^()]*\))?|\[([^\[\]]+(?=[\[\]])|\[[^\[\]]*\])*\])~/trans($1)/g;
    if ($in =~ /~/) {
      print "# Skipping (transpose notation) `$in'\nok $c\n";
      return;
    }
    if ($in !~ /\w\(/) {	# No function calls?
      # XXXX Primitive!
      # Substitute constants where they are not arguments to functions
      $in =~ s/(^|\G|\W)([-+]?\d+(\.\d*)?)(.?)/$1 PARI($2) $4/g;
      # Big integer constants:
      $in =~ s/\bPARI\((\d{10,})\)/PARI('$1')/g;
    } else {
      # Substitute big integer constants
      $in =~ s/(^|\G|\W)(\d{10,}(?!\.\d*))(.?)/$1 PARI('$2') $3/g;
      # Substitute division 
	$in =~ s/(^|[\-\(,\[])(\d+)\s*\/\s*(\d+)(?=$|[\),\]])/$1 PARI($2)\/PARI($3) /g;
    }
    # Substitute i= in loop commands
    if ($in !~ /\bhermite\(/) {	# Special case, not loop
      $in =~ s/([\(,]\w+)=(?!=)/$1,/g;
    }
    # Substitute print
    $in =~ s/\b(|p|tex)print\(/ 'my_' . $1 . 'print(1,' /ge;
    $in =~ s/\b(|p|tex)print1\(/ 'my_' . $1 . 'print(0,'/ge;
    $in =~ s/\b(eval|shift|sort)\(/&$1\(/g; # eval($y)
    # Recognize variables
    %seen_now = ();
    $in =~ s/(^|[;(])(\w+)(\s*=\s*)/$seen_now{$2} = '$'; $1 . '$' . $2 . $3/ge; # Assignment
    # Substitute variables (not before '^' - inside of 'o(x^17)'):
    $in =~ s/(^|[^\$])\b([a-zA-Z]\w*)\b(?!\s*[(^])/($1 || '') . ($seen{$2} || $seen_now{$2} || '') . $2/ge;
    # Die if did not substitute variables:
    while ($in =~ /(^|[^\$])\b([a-zA-Z]\w*)\b(?!\s*[\{\(^])/g) {
      print("# Skipping ($2 was not set) `$in'\nok $c\n"), return
	unless $seen{$2} and $seen{$2} eq ' ' or $in =~ /\"/;
      # Let us hope that lines which contain '"' do not contain unset vars
    }
    # Simplify for the following conversion:
    $in =~ s/\brandom\(\)/random/g;
    # Sub-ify sum,prod
    1 while
      $in =~ s/
		(
		  \b 
		  (
		    solve 
		  |
		    (?:
		      post (?! plothraw \b)
		    )?
		    ploth \w* 
		  |
		    sum \w* 
		  |
		    prod \w* 
		  |
		    v? vector 
		  |
		    matrix 
		  |
		    intgen 
		  |
		    intnum 
		  |
		    intopen 
		  |
		    for \w+
		  )
		  \( 
		  (?:
		    [^()]+ 
		    (?=
		      [(,)]
		    )
		  |
		    \( [^()]+ \)
		  )+		# One level of parenths supported
		  ,
		)
		(?!\s*sub\s*\{)	# Skip already converted...
		(		# This follows after a comma on toplevel
		  (?:
		    [^(,)\[\]]+ 
		    (?=
		      [()\[\]]
		    )
		  |
		    \(		# One level of parenths
		    (?:
		      [^()]+ 
		      (?=
			[()]
		      )
		    |
		      \( [^()]+ \) # Second level of parenths
		    )*
		    \)
		  |
		    \[		# One level of brackets
		    (?:
		      [^\[\]]+ 
		      (?=
			[\[\]]
		      )
		    |
		      \[ [^\[\]]+ \] # Second level of brackets
		    )*
		    \]
		  )*		# Two levels of parenths supported
		)
		\)
	      /$1 sub{$3}\)/xg;
    # Convert 10*20 to integer
    $in =~ s/(\d+)(?=\*\*\d)/ PARI($1) /g;
    # Convert blah[3], blah()[3] to blah->[-1+3]
    $in =~ s/([\w\)])\[/$1 -> [-1+/g;
    # Workaround for &eval test:
    $in =~ s/\$y=\$x;&eval\b(.*)/PARI('y=x');&eval$1;\$y=\$x/;
    # Workaround for kill:
    $in =~ s/^kill\(\$(\w+)\);/kill('$1');\$$1=PARIvar '$1';/;
    print "# eval", ($noans ? "-noans" : '') ,": $in\n";
    $printout = '';
    my $have_floats = ($in =~ /\d+\.\d*|\d{10,}/ 
		       or $in =~ /\b(zeta|bin|comprealraw|frac|lseriesell|powrealraw|legendre|suminf|forstep)\b/);
    # Remove the value from texprint:
    pop @$out if $in =~ /texprint/ and @$out == 2;
    $res = eval "$in";
    $rres = $res;
    $rres = pari_print $res if defined $res and ref $res;
    if ($doprint) {
      $rout = join "\t", @$out, "";
      if ($have_floats) {
	$printout = massage_floats $printout, "14f";
	$rout = massage_floats $rout, "14f";
      }
    } else {
      $rout = mformat @$out;
      if (not $doprint and $rout =~ /\[.*[-+,]\s/) {
	$rout =~ s/,* +/ /g;
	$rres =~ s/,* +/ /g if defined $res;
      }
    }
    if ($have_floats and ref $res) {
      if ($in =~ /\b(zeta|bin|comprealraw|frac|lseriesell|powrealraw|legendre|suminf)\b/) {
	$rres = massage_floats $rres, "14f";
	$rout = massage_floats $rout, "14f";
      } else {
	$rres = massage_floats $rres;
	$rout = massage_floats $rout;
      }
    }
    if ($@) {
      print "not ok $c # in='$in', err='$@'\n";
    } elsif (not $noans and (not defined $rres or $rres ne $rout)) {
      print "not ok $c # in='$in'\n#    out='", $rres, 
      "'\n# expect='$rout', type='", ref $res,"'\n";
    } elsif ($doprint and $printout ne $rout) {
      print "not ok $c # in='$in'\n# printout='", $printout, 
      "'\n#   expect='$rout', type='", ref $res,"'\n";
    } else {
      print "ok $c\n";
      @seen{keys %seen_now} = values %seen_now;
    }
  }
}

sub process_error {
  my ($in, $out, $error) = @_;
  $c++;
  print("# Skipping error($error) test $c: `$in'\nok $c\n");
}

sub process_definition {
  my ($name, $def) = @_;
  $c++;
  eval "PARI('$def');  import Math::Pari $name;";
  if ($@) {
    chomp $@;
    print("not ok $c # definition: `$in' error `$@'\n");
  } else {
    print("# definition $c: `$in'\nok $c\n");
  }
}

sub process_set {
  my ($in, $out) = @_;
  return process_test("setprecision($1)", 1, '') if $in =~ /^\\precision\s*=\s*(\d+)$/;
  $c++;
  print("# Skipping setting test $c: `$in'\nok $c\n");
}

sub process_print {
  my ($in, @out) = @_;
  $c++;
  print("# Skipping print $c: `$in'\nok $c\n");
}

sub process_multi {
  my ($in, $out) = @_;
  my @out = @$out;
  $c++;
  print("# Skipping multiline $c: `$in'\nok $c\n");
}

sub my_print {
  my $nl = shift;
  @_ = map {(ref) ? (pari_print $_) : $_} @_;
  $printout .= join '', @_;
  $printout .= "\t" if $nl;
  return;
}

sub my_pprint {
  my $nl = shift;
  @_ = map {(ref) ? (pari_pprint $_) : $_} @_;
  $printout .= join '', @_;
  $printout .= "\t" if $nl;
  return;
}

sub my_texprint {
  my $nl = shift;
  @_ = map {(ref) ? (pari_texprint $_) : $_} @_;
  $printout .= join '', @_;
  $printout .= "\t" if $nl;
  return;
}

