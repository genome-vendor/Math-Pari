# $Id: Makefile.PL,v 1.1 1995/01/23 04:46:12 ilya Exp $
use ExtUtils::MakeMaker;
use Config '%Config';

$paridir = "./libPARI";		# Change this line if you have it
$os = $Config{osname};

sub process_sparc {
  my $info = shift;
  return 'sparcv8_super' if $info =~ /TMS390Z5[05]|CY605/;
  return 'sparcv8_micro' if $info =~ /TMS390S10|MB8690[47]/;
  return 'sparcv8_hyper' if $info =~ /RT625/;
  return shift;
}

sub process_machine {
  my $os = (split ' ', $Config{myuname})[0];
  my $machine = $os;
  if ($os =~ /^hp/) {
    $machine = 'hppa';
  } elsif ($os eq 'os2' or $os eq 'freebsd') {
    $machine = 'i386';
  } elsif ($os eq 'ultrix') {
    $machine = 'mips';
  } elsif ($os eq 'nextstep' or -d '/NextApps') {
    chomp($machine = `file /bin/sh | sed 's/.*(for architecture \(.*\))/\1/'`);
  } elsif ($os eq 'osf1') {
    $machine = 'alpha' if (split ' ', $Config{myuname})[4] eq 'alpha';
  } elsif ($os eq 'linux') {
    chomp($machine = `uname -m`);
    if (-e '/proc/cpuinfo') {
      open IN, '/proc/cpuinfo' or die "open /proc/cpuinfo: $!";
      local $/ = undef;		# Needed?
      my $info = <IN>;
      close IN or die "close /proc/cpuinfo: $!";
      $machine = process_sparc $info, $machine;
    }
  } elsif ($os eq 'sunos') {
    my $type = (split ' ', $Config{myuname})[4];
    if ($type =~ /^sun3/) {
      $machine = 'm68k';
    } elsif ($type =~ /^sun4[ce]/) {
      $machine = 'sparcv7';
    } elsif ($type =~ /^sun4[dm]/) {
      local $ENV{PATH} = "$ENV{PATH}:/dev/sbin";
      my $info = `(prtconf||devinfo)2>&-`;
      $info = join ' ', grep /TI,|FMI,|Cypress,|Ross,/, split "\n", $info;
      $machine = process_sparc $info, $machine;
    } elsif ($type eq 'sun4u') {
      $machine = 'sparcv9';
    }
  }
  print "Processor $machine detected\n";
  # For older PARI:
  $machine = 'sparcv8super'
    if $machine eq 'sparcv9' or $machine eq 'sparcv8_hyper'
       or $machine eq 'sparcv8_super';
  $machine = 'sparcv8micro' if $machine eq 'sparcv8_micro';
  return $machine;
}

if (@ARGV and $ARGV[0] =~ /^machine=(.*)$/i) {
  $common::machine = $1; shift;
} else {
  $common::machine = process_machine;
}

$define = "-Dshifts=pari_shifts -DDEBUG_PARI";
#$libs = "-L$paridir -lpari -lgnuplot -lm";
$libs = "-lgnuplot -lm";

# Comment if your Pari is not c++, or your compiler is not gcc:
$cpp = "-x c++";

if ($os eq "solaris") {
  # Below -R specifies the path to search for this file when the program will run, -L the path to search when the program is linked
  @sc_dirs = </opt/SUNWspro/SC*>;
  die "Cannot find SUNWspro dirs, needed for -lsunmath" unless @sc_dirs;
  $libs .= " -L$sc_dirs[-1]/lib -lsunmath -lm";
} elsif ($os eq "os2") {
  $define = "$define -DNOEXP2 -DMALLOC_PROCS";
} elsif ($os eq 'linux') {
  $define .= ' -DNOEXP2';
}


# We want to find the offset of the least significant byte
# of SV.flags inside an SV.
# It is enough to find the offset of the least significant byte inside a U32.
# (Needed for recognizing Perl code supplied to lisseq and lisexpr.)

$offset = index($Config{byteorder}, "1");
$offset = 3 if $offset > 3;	# SvFLAGS is actually I32

$define .= " -DLSB_in_U32=$offset";

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile being created.
&WriteMakefile(
	       'LIBS' =>  $libs,
	       'INC' => "-I$paridir",
	       'NAME' => 'Math::Pari',
	       'MYEXTLIB' => 'libPARI/libPARI$(LIB_EXT)',
	       'VERSION_FROM' => 'Pari.pm',
	       'DEFINE' => $define,
	       'dist' => {COMPRESS=>'gzip -9f', SUFFIX=>'gz'},
	       'clean'      => { FILES => 'pari.ps' },
	      );


sub MY::postamble {
  '
$(MYEXTLIB): libPARI/Makefile FORCE
	cd libPARI && $(MAKE) $(PASTHRU)
';
}
 
