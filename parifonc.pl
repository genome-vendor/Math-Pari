unshift @ARGV, 'Pari.xs' if @ARGV < 2 and -r 'Pari.xs';
unshift @ARGV, 'libPARI/anal.c' if @ARGV < 2 and -r 'libPARI/anal.c';
(@ARGV == 2) || &usage;

@known = split /,\s*/, 'label, while, goto, until, read, pprint, print, texprint, pprint1, print1, O, if, o';
@known{@known} = (1) x @known;

open(XSUB,$ARGV[1]) || die "Cannot open $ARGV[1]: $!";
while (<XSUB>) {
  $supported{$1}++ if /^\s*CASE_INTERFACE\s*\(\s*(\d+)\s*\)/;
}
close(XSUB) || die "Cannot close $ARGV[1]: $!";

open(ANAL,$ARGV[0]) || die "Cannot open $ARGV[0]: $!";
while (<ANAL>) {
  if (/^entree\s+fonctions\[/ ... /^\s*\}\s*;\s*$/) {
    next unless $i++;		# Skip first line
    last if /^\s*\}\s*;\s*$/;
    &warnl() unless /
		      ^ \s* \{ \s* " 
		      (
			[^""]+	# 1 Name
		      )
		      " \s* , \s* 
		      (
			\d+	# 2 Interface
		      )
		      \s* , \s* 
		      (
			[^,]+	# 3 C function pointer
		      )
		      \s* , \s* 
		      (
			\d+	# 4 Group
		      )
		      \s* , \s* 
		      (
			\d+	# 5 
		      )
		      ( \s* , \s* ((" [^"]* ") | NULL) \s* , \s* NULL )? # New fields
		      \s* \} \s* ,? \s* $
		    /x; # ";
    ($pari, $interface, $gp, $group, $code) = ($1, $2, $3, $4, $8);
    if ($gp eq "0") {
      if ($known{$pari}) {
	$builtin_known{$pari}++;
      } else {
	$builtin{$pari}++;
      }
    } else {
      $interface{$pari} = $interface;
      $code{$interface} ||= ($code || '');
      push @{$group{$group}}, $pari unless exists $supported{$interface};
      $interfaces{$interface}++;
    }
    # print "'$pari' <= '$gp' via $interface\n";
    # &warnl() unless $gp =~ /\b$pari\b/;
  }
}
close(ANAL) || die "Cannot close $ARGV[0]: $!";
print "Builtins, unsupported as functions (but available in Perl):\n\t", join(", ", keys %builtin_known), "\n\n"
  if %builtin_known;

print "Builtins, completely unsupported:\n\t", join(", ", keys %builtin), "\n\n"
  if %builtin;

for (keys %interfaces) {
  $unsupported{$_}++ unless $supported{$_};
}

@unsupported = sort {$interfaces{$a} <=> $interfaces{$b}} keys %unsupported;

print "\tTotal number of unsupported interfaces: ",scalar @unsupported,":\n";
for $i (sort {$a <=> $b} @unsupported) {
  print "Interface $i$code{$i} used in $interfaces{$i} function(s): ",
     join(", ", @f=grep($interface{$_}==$i, keys %interface)), ".\n";
  $total += $interfaces{$i};
  push(@ff,@f);
}

print "\n\tTotal number of unsupported functions: $total:\n";
  #join(", ", sort @ff), "\n";

for $g (sort {$a <=> $b} keys %group) {
  print "group $g:\t", join(', ', sort @{$group{$g}}), "\n";
}

sub usage {die "Usage: $0 [path/to/anal.c] [path/to/Pari.xs]\n"}

sub warnl {warn "Unrecognized line:\n$_"}

