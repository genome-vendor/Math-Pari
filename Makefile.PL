# $Id: Makefile.PL,v 1.1 1995/01/23 04:46:12 ilya Exp $
use ExtUtils::MakeMaker;
use Config '%Config';

$paridir = "./libPARI";		# Change this line if you have it
$os = $Config{osname};

if (@ARGV and $ARGV[0] =~ /^machine=(.*)$/i) {
  $common::machine = $1; shift;
}

$define = "-DDEBUG_PARI";
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
	       'VERSION' => '0.8',
	       'DEFINE' => $define,
	       'dist' => {COMPRESS=>'gzip -9f', SUFFIX=>'gz'},
	      );


sub MY::postamble {
  '
$(MYEXTLIB): libPARI/Makefile
	cd libPARI && $(MAKE) $(PASTHRU)
';
}
 
