# $Id: Makefile.PL,v 1.1 1995/01/23 04:46:12 ilya Exp $
use ExtUtils::MakeMaker;
use Config '%Config';
use File::Copy 'copy';

# Uncomment and change this line if I fail to find it
# $paridir = "..";

$latmus = 'src/test/in/nfields';

sub format_pari {
  my $dir = shift;
  $dir =~ s/(\d+)\.(\d+).(\d+)/sprintf('%03d%03d%03d',$1,$2,$3)/e;
  $dir;
}

unless (defined $paridir) {
  # Assume we are inside:
  for $dir ('..', '../..') {
    $paridir = $dir, last 
      if -e "$dir/$latmus";
  }
  unless (defined $paridir) {
    # Try to find alongside
    for $dir ('.', '..', '../..') {
      @dirs = <$dir/pari-[234].*>;
      @dirs = grep -e "$_/$latmus", @dirs;
      last if @dirs;
    }
    @gooddirs = grep !/alpha|beta/, @dirs;
    @gooddirs = grep !/alpha/, @dirs unless @gooddirs;
    @gooddirs = @dirs unless @gooddirs;
    @gooddirs = sort {format_pari($a) cmp format_pari($b)} @gooddirs;
    $paridir = $gooddirs[-1] if @gooddirs;
  }
}
if ($paridir) {
  print <<EOP;
Found PARI build directory in $paridir 
     (formatted @{[format_pari($paridir)]}).
EOP
} else {
  die "Could not find PARI build directory, please edit Makefile.PL.\n";
}
$common::main_paridir = $paridir;
$os = $Config{osname};

if (@ARGV and $ARGV[0] =~ /^machine=(.*)$/i) {
  $common::machine = $1; shift;
}

$define = "-Dshifts=pari_shifts -DDEBUG_PARI";
#$libs = "-L$paridir -lpari -lgnuplot -lm";
$libs = "-lgnuplot -lm";

# Comment if your Pari is not c++, or your compiler is not gcc:
# $cpp = "-x c++";

$noexp2 = '';			# Should be used for generation of paricfg.h

if ($os eq "solaris") {
  # Below -R specifies the path to search for this file when the
  # program will run, -L the path to search when the program is linked
  @sc_dirs = </opt/SUNWspro/SC*>;
  if (@sc_dirs) {
    $libs .= " -L$sc_dirs[-1]/lib -lsunmath -lm";
  } else {
    warn "Cannot find SUNWspro dirs, needed for -lsunmath, using NOEXP2.\n";
    $noexp2 = 1;
  }
} elsif ($os eq "os2") {
  $noexp2 = 1;
  $define .= " -DMALLOC_PROCS";
} elsif ($os eq 'linux') {
  $noexp2 = 1;
}


# We want to find the offset of the least significant byte
# of SV.flags inside an SV.
# It is enough to find the offset of the least significant byte inside a U32.
# (Needed for recognizing Perl code supplied to lisseq and lisexpr.)

$offset = index($Config{byteorder}, "1");
$offset = 3 if $offset > 3;	# SvFLAGS is actually I32

$define .= " -DLSB_in_U32=$offset";
$parisrc = "$paridir/src";
$paritests = "$parisrc/test/in";

opendir TESTS, $paritests 
  or die "Cannot find tests in $paritests: $!";
@tests = readdir TESTS;
closedir TESTS or die "Cannot find tests (close): $!";
for $test (@tests) {
  next if $test =~ /^\.\.?$/;
  next if $test =~ /compat/;
  $targ = "t/$test.t";
  $sou = 'test_eng/ex.t';
  next if -e $targ and -M $targ <= -M $sou;
  if ($paridir ne '..' and not -e "$sou-") {
    system "$^X -pe 's,\.\.,$paridir,' $sou > $sou- "
      and die "Could not run test converter: $! $?";
    $sou = "$sou-";
  }
  if (-f $targ) {
    chmod 0666, $targ;
    unlink $targ;
  }
  copy $sou, $targ or die "Cannot create test $test.t: $1";
}
$targ = 'libPARI.pod';
if (not -e $targ
    or -M $targ > -M "$paridir/doc/usersch3.tex"
    or -M $targ > -M "paridoc_to_pod") {
  if (-f $targ) {
    chmod 0666, $targ;
    unlink $targ;
  }
  system "$^X paridoc_to_pod $paridir/doc/usersch3.tex > $targ "
    and die "Errors when converting documentation: $! $?"
}

@paricfg = <$paridir/o.*/paricfg.h>;
if (@paricfg > 1) {
  warn "Found multiple paricfg.h: @paricfg.\n";
  @paricfg = sort { -M $a <=> -M $b} @paricfg;
  warn "Choosing newest paricfg.h: $paricfg[0].\n";
} 
if (@paricfg == 0) {
  warn <<EOW;
Did not find paricfg.h.  You will need to manually copy it to libPARI
    directory from the PARI build directory.
EOW
} else {
  copy $paricfg[0], 'libPARI/paricfg.h' 
    or die "Could not copy $paricfg[0] to paricfg.h: $!"
      if not -e 'libPARI/paricfg.h' or -M $paricfg[0] < -M 'libPARI/paricfg.h';
}

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile being created.
&WriteMakefile(
	       LIBS	     =>  $libs,
	       INC	     => "-I$parisrc/headers -I$parisrc -I./libPARI",
	       NAME	     => 'Math::Pari',
	       MYEXTLIB	     => 'libPARI/libPARI$(LIB_EXT)',
	       VERSION_FROM  => 'Pari.pm',
	       DEFINE	     => $define,
	       dist	     => {COMPRESS=>'gzip -9f', SUFFIX=>'gz'},
	       clean	     => { FILES => 'pari.ps' },
	      );


sub MY::postamble {
  '
$(MYEXTLIB): libPARI/Makefile FORCE
	cd libPARI && $(MAKE) $(PASTHRU)
';
}
 
