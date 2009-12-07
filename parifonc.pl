(@ARGV == 2) || &usage;
open(ANAL,$ARGV[0]) || die "Cannot open $ARGV[0]: $!";
while (<ANAL>) {
  if (/^entree\s+fonctions\[/ ... /^\s*\}\s*;\s*$/) {
    next unless $i++;		# Skip first line
    last if /^\s*\}\s*;\s*$/;
    &warnl() unless /^\s*\{\s*"([^""]+)"\s*,\s*(\d+)\s*,\s*([^,]+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\}\s*,?\s*$/; # ";
    ($pari, $interface, $gp) = ($1, $2, $3);
    if ($gp eq "0") {
      $builtin{$pari}++;
    } else {
      $interface{$pari} = $interface;
      $count{$interface}++;
    }
    # print "'$pari' <= '$gp' via $interface\n";
    # &warnl() unless $gp =~ /\b$pari\b/;
  }
}
close(ANAL) || die "Cannot close $ARGV[0]: $!";
print "Builtins, unsupported:\n\t", join(", ", keys %builtin), "\n\n";

open(XSUB,$ARGV[1]) || die "Cannot open $ARGV[1]: $!";
while (<XSUB>) {
  $supported{$1}++ if /^\s*CASE_INTERFACE\s*\(\s*(\d+)\s*\)/;
}
close(XSUB) || die "Cannot close $ARGV[1]: $!";

for (keys %count) {
  $unsupported{$_}++ unless $supported{$_};
}

@unsupported = sort {$count{$a} <=> $count{$b}} keys %unsupported;

for $i (@unsupported) {
  print "Unsupported interface $i used in $count{$i} function(s): ",
     join(", ", @f=grep($interface{$_}==$i, keys %interface)), ".\n";
  $total += $count{$i};
  push(@ff,@f);
}

print "\n\tTotal number of unsupported functions: $total:\n",
  join(", ", sort @ff), "\n";

sub usage {die "Usage: $0 path/to/anal.c path/to/Pari.xs\n"}

sub warnl {warn "Unrecognized line:\n$_"}

