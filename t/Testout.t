#! perl
# 	$rcs = ' $Id: testout.t,v 1.2 1997/09/22 10:13:37 ilya Exp ilya $ ' ;	

use Math::Pari qw(:DEFAULT pari_print :all);
use vars qw($x $y $z $k $t $q $a $u $j $l $name $other);
print "1..565\n";

push @ARGV, 'libPARI/testouts' unless @ARGV;
$| = 1;
@seen{qw(pi i euler a x y z k t q u j l name other)} = (' ', ' ', ' ', ('$') x 100);
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
  $wasprint = /\b((|p|tex)print|plot)\b/;
  chomp;
  $in = $_;
  $_ = '';
  $_ = <> while defined and /^$/; # Before warnings
  $_ = <> while defined and /^\s*\*+\s*warning/i; # skip warnings
  if (s/^\s*\*+\s*//) {		# error
    process_error($in,$_);
    next;
  }
  defined or die "Can't find an answer";
  chomp;
  process_test($in, 1, ''), redo if /^\?\s/; # Was a void 
  s/^%\d+\s*=\s*// or die "Malformed answer: $_" unless $bad or $wasprint;
  if ($_ eq '' or $wasprint) {	# Answer is multiline
    @ans = ();
    while (<>) {
      last if /^\?\s+/;
      next if /^$/;
      chomp;
      push @ans, $_;
    }
    if ($wasprint) {
      process_print($in, @ans);
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

sub process_test {
  my ($in, $noans, $out) = @_;
  $c++;
  # First a trivial processing:
  $in =~ s/\b(\d+)\s*\\\s*(\d+)/ gdivent($1,$2)/g; # \
  $in =~ s/\b(\d+)\s*\\\/\s*(\d+)/ gdivround($1,$2)/g;	# \/
  $in =~ s/\b(\d+)\s*!/ ifact($1)/g; # !
  if ($in =~ /\\/) {		# \\ for division unsupported
    $c--;
    process_error($in, $out);
  } elsif ($in =~ /^(\w+)\s*\([^()]*\)\s*=/) {
    $c--;
    process_definition($1, $in);    
  } elsif ($in =~ /!/) { # Factorial
      print "# Skipping (factorial) `$in'\nok $c\n";
  } else {
    # work with "^", need to treat differently inside o()
    $in =~ s/\^/^^^/g;
    $in =~ s/\bo\(([^()]*)\^\^\^([^()]*)\)/ PARI('O($1^$2)') /gi;
    $in =~ s/\^\^\^/**/g;	# Now treat it outside of O()
    $in =~ s/\[([^\[\];]*;[^\[\]]*)\]/format_matrix($1)/ge; # Matrix
    $in =~ s/\[([^\[\];]*)\]\s*~/format_vvector($1)/ge; # Vertical vector
    if ($in =~ /\[[^\]]*;/) { # Matrix
      print "# Skipping (matrix notation) `$in'\nok $c\n";
      return;
    } elsif ($in =~ /(^|[\(=])%/) {
      print "# Skipping (history notation) `$in'\nok $c\n";
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
      $in =~ s/(^|[\(,])(\d+)\s*\/\s*(\d+)($|[\),])/$1 PARI($2)\/PARI($3) $4/g;
    }
    if ($in =~ /\bget(heap|stack)/) { # Meaningless
      print "# Skipping meaningless `$in'\nok $c\n";
    } elsif ($in =~ /\bsettype\b/) {
      # Dies XXXXX
      print "# Skipping (possibly FATAL $1) `$in'\nok $c\n";
    } elsif ($in =~ /\b(v?vector\b|matrix\b|settype\b)/) {
      # Dies on settype XXXXX
      print "# Skipping (converting $1 needs additional work) `$in'\nok $c\n";
    } else {
      # Recognize variables
      %seen_now = ();
      $in =~ s/(^|[;(])(\w+)(\s*=\s*)/$seen_now{$2} = '$'; $1 . '$' . $2 . $3/ge; # Assignment
      # Substitute variables (not before '^' - inside of 'o(x^17)'):
      $in =~ s/(^|[^\$])\b([a-zA-Z]\w*)\b(?!\s*[(^])/($1 || '') . ($seen{$2} || $seen_now{$2} || '') . $2/ge;
      # Die if did not substitute variables:
      while ($in =~ /(^|[^\$])\b([a-zA-Z]\w*)\b(?!\s*[(^])/g) {
	print("# Skipping ($2 was not set) `$in'\nok $c\n"), return
	  unless $seen{$2} and $seen{$2} eq ' ';
      }
      print "# eval", ($noans ? "-noans" : '') ,": $in\n";
      $res = eval "$in";
      $out = mformat @$out;
      if ($@) {
	print "not ok $c # in='$in', err='$@'\n";
      } elsif (not $noans and ref $res and pari_print($res) ne $out) {
	print "not ok $c # in='$in', out='", pari_print($res), 
	"', expect='$out', type='", ref $res,"'\n";
      } elsif (not $noans and not ref $res and $res ne $out) {
	print "not ok $c # in='$in', out='", pari_print($res), 
	"', expect='$out', type='", ref $res,"'\n";
      } else {
	print "ok $c\n";
	@seen{keys %seen_now} = values %seen_now;
      }
    }
  }
}

sub process_error {
  my ($in, $out) = @_;
  $c++;
  print("# Skipping error test $c: `$in'\nok $c\n");
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

