# $Id: Makefile.PL,v 1.1 1995/01/23 04:46:12 ilya Exp $
use ExtUtils::MakeMaker;
use Config '%Config';
use File::Copy 'copy';

# Uncomment and change this line if I fail to find it
# $paridir = "..";

unless (-t STDOUT) {		# Better reports if run with redirections
    my $odef = select STDERR;
    $| = 1;
    select STDOUT;
    $| = 1;
    select $odef;  
}

$latmus = 'src/test/in/nfields';

sub get_pari_version {
  my $dir = shift;
  my $v = "";
  open(IN, "$dir/config/version") or return;
  /(?:VersionMajor|VersionMinor|patch)='?(\d+(\.\d+)?)'?/
    and $v .= "$1." while <IN>;
  close(IN) or die "error closing '$dir/config/version'";
  $v =~ s/\.$// or return;
  return "pari-$v";
}

sub format_pari {
  my $dir = shift;
  my $v;
  $v = (get_pari_version $dir
	or warn(<<EOW), $dir);
Could not extract version from '$dir/config/version', using directory name...
EOW
  $dir =~ s,/pari$,/$v,;	# For CVS 'pari' directory
  $dir =~ s/(\d+)\.(\d+).(\d+)/sprintf('%03d%03d%03d',$1,$2,$3)/e;
  $dir;
}

unless (defined $paridir) {
  # Assume we are inside:
#  for $dir ('..', '../..') {
#    $paridir = $dir, last 
#      if -e "$dir/$latmus";
#  }
  unless (defined $paridir) {
    # Try to find alongside
    for $dir ('.', '..', '../..', '../../..') {
      @dirs = <$dir/pari-[234].*>;
      @dirs = "$dir/pari" if not @dir and -d "$dir/pari";
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

$paridir = get_pari() if !$paridir;  

if ($paridir) {
  $pari_version = format_pari($paridir);
  print <<EOP;
Found GP/PARI build directory in $paridir 
     (formatted $pari_version).
EOP
  $pari_version =~ s/.*pari-0+([0-9]*).*/$1/;
  $common::pari_version = $pari_version;
} else {
  die "Could not find GP/PARI build directory, please edit Makefile.PL.\n";
}
$common::main_paridir = $paridir;
$os = $Config{osname};

grep {/^machine=(.*)$/i and $common::machine = $1} @ARGV;
grep {/^configure$/i and $do_configure = 1} @ARGV;

@ARGV = grep !/^(machine=|configure$)/i, @ARGV;

#$define = "-Dshifts=pari_shifts -DDEBUG_PARI";
$define = "-DPARI_VERSION_EXP=$pari_version -DDEBUG_PARI";
#$libs = "-L$paridir -lpari -lgnuplot -lm";
#$libs = "-lgnuplot -lm";
$libs = "-lm";

# Comment if your Pari is not c++, or your compiler is not gcc:
# $cpp = "-x c++";

$noexp2 = '';			# Should be used for generation of paricfg.h

if ($os eq "solaris") {
  # Below -R specifies the path to search for this file when the
  # program will run, -L the path to search when the program is linked
  @sc_dirs = '/opt/SUNWspro/lib' if -d '/opt/SUNWspro/lib';
  @sc_dirs = </opt/SUNWspro/SC*/lib> unless @sc_dirs;
  if (@sc_dirs) {
    $libs .= " -R$sc_dirs[-1] -L$sc_dirs[-1] -lsunmath -lm";
  } else {
    #warn "Cannot find SUNWspro dirs, needed for -lsunmath, using NOEXP2.\n";
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
$sou = 'test_eng/ex.t';
next if -e $targ and -M $targ <= -M $sou;
#if (not -e "$sou-" or -M "$sou-" > -M $sou) {
  system "$^X -pe 's,CHANGE_ME,$paridir,' $sou > $sou- "
    and die "Could not run test converter: $! $?";
#}
$sou = "$sou-";
for $test (@tests) {
  next if $test =~ /^\.\.?$/;
  next if $test =~ /compat/;
  next if -d "$paritests/$test" and $test eq 'CVS';
  next if $test =~ /(~|\.(bak|orig|rej))$/;
  $targ = "t/$test.t";
  if (-f $targ) {
    chmod 0666, $targ;
    unlink $targ;
  }
  copy $sou, $targ or die "Cannot create test $test.t: $1";
}

# We can't do what the commented chunk does: $paridir/doc/gphelp is
# auto-generated, so its date is not relevant to anything.

# $targ = 'libPARI/gphelp';
# if (not -e $targ
#     or -M $targ > -M "$paridir/doc/gphelp") {
#   if (-f $targ) {
#     chmod 0666, $targ;
#     unlink $targ;
#   }
#   copy "$paridir/doc/gphelp", $targ;
# }

sub make_pod {
  my ($targ, $how) = @_;
  if (not -e $targ
      or -M $targ > -M "$paridir/doc/usersch3.tex"
      or -M $targ > -M "libPARI/gphelp") {
    if (-f $targ) {
      chmod 0666, $targ;
      unlink $targ;
    }
    (system "$^X libPARI/gphelp $how $paridir/doc/usersch3.tex > tmp_pod "
      and (warn("Errors when converting documentation: $?"), 0))
      or rename 'tmp_pod', $targ;
  }
}

make_pod 'libPARI.pod', '-to_pod';
make_pod 'libPARI.dumb.pod', '-to_dumb_pod';

@paricfg = <$paridir/o.*/paricfg.h>;
push @paricfg, <$paridir/O*/paricfg.h>;

# Reported to work with Win32 too
@paricfg = grep !/Odos/, @paricfg unless $^O =~ /dos|djgcc|MSWin32/i;

# Probably not present in newer versions anymore...
unshift @paricfg, "$paridir/win32/paricfg.h"
  if $^O eq 'MSWin32' and -f "$paridir/win32/paricfg.h";

if (not @paricfg and $do_configure) {
  print  "cd $paridir ; sh ./Configure\n";
  system "cd $paridir ; sh ./Configure"
    and die "Cannot configure: $!, exitcode=$?.\n";
  print "Configuration of GP/PARI successful.\n";
  @paricfg = <$paridir/O*/paricfg.h>;

  # Reported to work with Win32 too
  @paricfg = grep !/Odos/, @paricfg unless $^O =~ /dos|djgcc|MSWin32/i;
}

if ($do_configure) {
  if (@paricfg > 1) {
    warn "Found multiple paricfg.h: @paricfg.\n";
    @paricfg = sort { -M $a <=> -M $b} @paricfg;
    warn "Choosing newest paricfg.h: $paricfg[0].\n";
  }
  if (@paricfg == 0) {
    warn <<EOW;
Did not find paricfg.h.  You may need to manually copy it to libPARI
    directory from the GP/PARI build directory.
    ...Now switching to creation of paricfg.h by Perl code.
EOW
    $do_configure = 0;
  } else {
    copy $paricfg[0], 'libPARI/paricfg.h' 
      or die "Could not copy $paricfg[0] to paricfg.h: $!"
	if not -e 'libPARI/paricfg.h' or -M $paricfg[0] < -M 'libPARI/paricfg.h';
  }
}

my $have_ladd;
write_configure() unless -r 'libPARI/paricfg.h';

$define .= ' -Derr=pari_err'; # On linux it can get a wrong dynamic loading

$define .= ' -DHAVE_LADD' if $have_ladd;

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile being created.
&WriteMakefile(
	       LIBS	     =>  $libs,
	       INC	     => '-I$(PARI_DIR)/src/headers -I$(PARI_DIR)/src -I./libPARI',
	       NAME	     => 'Math::Pari',
	       MYEXTLIB	     => 'libPARI/libPARI$(LIB_EXT)',
	       VERSION_FROM  => 'Pari.pm',
	       DEFINE	     => $define,
	       dist	     => {COMPRESS=>'gzip -9f', SUFFIX=>'gz'},
	       clean	     => { FILES => 'pari.ps' },
	       macro	     => { 
				 PARI_DIR	  => $paridir,
				},
	      );


sub MY::postamble {
  '
$(MYEXTLIB): libPARI/Makefile FORCE
	cd libPARI && $(MAKE) $(PASTHRU)

#libPARI/gphelp: $(PARI_DIR)/doc/gphelp
#	$(CP) -f $(PARI_DIR)/doc/gphelp $@

libPARI.pod: $(PARI_DIR)/doc/usersch3.tex libPARI/gphelp
	$(PERL) libPARI/gphelp -to_pod $(PARI_DIR)/doc/usersch3.tex > pod_tmp && mv pod_tmp $@

libPARI.dumb.pod: $(PARI_DIR)/doc/usersch3.tex libPARI/gphelp
	$(PERL) libPARI/gphelp -to_dumb_pod $(PARI_DIR)/doc/usersch3.tex > pod_tmp && mv pod_tmp $@
';
}

sub MY::const_config
{
 my $self = shift;
 my $flags = $self->{'CCCDLFLAGS'}; # Tmp var needed with Perl4 !
 $flags =~ s/(-[fK]?\s*)pic\b/${1}PIC/; 
 $self->{'CCCDLFLAGS'} = $flags;
 if ($^O eq 'MSWin32' && $Config{'ccflags'} =~ /-DPERL_OBJECT/)
  {
   $self->{'LDFLAGS'} =~ s/-(debug|pdb:\w+)\s+//g;
   $self->{'LDDLFLAGS'} =~ s/-(debug|pdb:\w+)\s+//g;
  }
 return $self->MM::const_config;
}

my $have_ulong;
my $clk_tck_def;
my $have_getrusage;

sub scan_headers {
  warn "Scanning header files...\n";
  my $cmd = "$Config{cpprun} $Config{cppflags} utils/inc.h 2>&1";
  open INC, "$cmd |" or warn("Error $! from: $cmd\n"), return;
  $clk_tck_def = 1;
  while (<INC>) {
#    print;
    
    $have_ulong = 1, warn "...ulong\n" if /\btypedef\b.*\bulong\s*;/;
    $clk_tck_def = 0, warn "...CLK_TCK not defined\n" if /y\s*=\s*CLK_TCK\b/;
    $have_getrusage = 1, warn "...getrusage\n" if /\bgetrusage\s*\(/;
    $have_ladd = 1, warn "...ladd\n" if /\bladd\b/;
  }
  close INC or warn "Note (probably harmless): Errors reading from pipe: '$!', exit=$?: $cmd\n"
}

sub write_configure {
  scan_headers() or $clk_tck_def = 0;
  warn "Creating libPARI/paricfg.h...\n";
  open F, '> libPARI/paricfg.h' or die "open 'libPARI/paricfg.h' for write: $!";
  print F <<EOP;
#define UNIX
#define GPDATADIR "/usr/local/lib/pari/galdata"
#define GPMISCDIR "/usr/local/lib/pari"

#define PARI_BYTE_ORDER    $Config{byteorder}
#define NOEXP2	/* Otherwise elliptic.t:11 rounds differetly, and fails */

EOP
  if ($have_getrusage) {
    print F <<EOP if $Config{d_times};
#define USE_GETRUSAGE 1

EOP
  } else {
    print F <<EOP if $Config{d_times};
#define USE_TIMES 1

EOP
    print F <<EOP if $Config{d_times} and $Config{i_time} and !$clk_tck_def;
/* Reported to be needed on some Linuxes: */
#include <time.h>

EOP
    print F <<EOP if not $Config{d_times} and $Config{d_ftime};
#define USE_FTIME 1

EOP
  }
  print F <<EOP if $Config{dlsrc} eq 'dl_dlopen.xs';
#define HAS_DLOPEN

EOP

  print F <<EOP unless $have_ulong;
#define ULONG_NOT_DEFINED

EOP

  close F or die "close 'libPARI/paricfg.h' for write: $!";
}

sub get_pari {
  $host = 'megrez.math.u-bordeaux.fr';
  $dir  = '/pub/pari/unix/';

  print "Did not find GP/PARI build directory around.\n";

  if (-t STDIN and (-t STDOUT or -p STDOUT)) { # Interactive
    $| = 1;
    $mess = <<EOP;
Do you want to me to fetch GP/PARI automatically?  
Make sure you have a large scrollback buffer to see the messages.
Fetch? (y/n, press Enter)
EOP
    chomp $mess;
    print "$mess ";
    $ans = <>;
    if ($ans !~ /y/i) {
      print "Well, as you wish...  I give up...\n";
      return;
    }
  } else {
    print "Non-interactive session, autofetching...\n"
  }

  print "Getting GP/PARI from ftp://$host$dir\n";

  eval { require Net::FTP }
    or die "You do not have Net::Ftp installed, cannot download, exiting...";

  $ftp = Net::FTP->new($host) or die "Cannot create FTP object: $!";
  $ftp->login("anonymous","Math::Pari@")
    or die "Cannot login anonymously (",$ftp->message(),"): $!";
  $ftp->cwd($dir) or die "Cannot cwd (",$ftp->message(),"): $!";
  $ftp->binary() or die "Cannot switch to binary (",$ftp->message(),"): $!";
  @lst = $ftp->ls() or die "Cannot list (",$ftp->message(),"): $!";
  #print "list = `@lst'\n";

  $match = 'pari.*(\d+\.\d+\.\d+).*\.t(ar\.)?gz$';
  $c = 0;
  for $file (@lst) {
    next unless $file =~ /$match/o;
    $c++;
    $version = $1;
    if ($file =~ /alpha/) {
      $alpha{$version} = $file;
    } elsif ($file =~ /beta/) { 
      $beta{$version} = $file;
    } else {
      $golden{$version} = $file;
    }
  }
  die "Did not find any file matching /$match/ via FTP"
    unless $c;

  sub fmt_version {sprintf "%03d%03d%03d", split /\./, shift}

  for $type (qw(alpha beta golden)) {
    if (%$type) {
      $have{$type}++;
      $best = $type;
      @files = keys %$type;
      print "Available ${type} versions: `@files'\n";
      $ {"latest_$type"} = (sort {fmt_version($a) cmp fmt_version($b)} 
			    keys %$type)[-1];
      $ {"latest_${type}_file"} = $$type{$ {"latest_$type"}};
      print qq(Latest $type is ${"latest_${type}_file"}\n);
    }
  }

  # Special-case v2.0.14
  if (!%golden and $latest_beta eq '2.0.11' and $latest_alpha eq '2.0.14') {
    $best = 'alpha';		# It is tested!
  }

  if ($best) {
    $file = $ {"latest_${best}_file"};
    $version = $ {"latest_$best"};
    print qq(Picking $best version $version, file $file\n);
    if (-f $file) {
      print qq(Well, I already have it, using the disk copy...\n);
    } else {
      print qq(Downloading...\n);
      $ftp->get($file) or die "Cannot get (",$ftp->message(),"): $!";
      print qq(Downloaded...\n);
    }
    print qq(Extracting...\n);
    my $zcat = "gzip -dc";	# zcat may be the old .Z-extractor
    print  "$zcat $file | tar -xvf -\n";
    system "$zcat $file | tar -xvf -"
      and die "Cannot extract: $!, exitcode=$?.\n";
    ($dir = $file) =~ s/\.t(ar\.)?gz$// or die "malformed name `$file'";
    -d $dir or die "Did not find directory $dir!";
    # Use GP/PARI's Configure to build paricfg.h
    unless ($^O =~ /dos|djgcc|MSWin32/i) { # We use an existing paricfg.h...
      if ($do_configure) {
	print  "cd $dir ; sh ./Configure\n";
	system "cd $dir ; sh ./Configure"
	  and die "Cannot configure: $!, exitcode=$?.\n";
	print "Configuration of GP/PARI successful.\n";
      } else {
	write_configure();
      }
    }
    my %patches = ('2.0.11'
		   => [qw(
			 patch11/diff_pari_gnuplot_aa
			 patch11/patch_pari_round0
			 patch11/patches_round1_short
			 patch11/diff_pari_fixed_interfaces_011
			 patch11/diff_pari_highlevel_hash_011a
			 patch11/diff_pari_ret_proto_2011)],
		   '2.0.12' => ['patch12/diff_for_perl_2012'],
		   '2.0.13' => ['patch13/diff_for_perl_2013',
				'patch13/diff_for_gnuplot_2013'],
		   '2.0.14' => ['patch14/diff_for_perl_2014',
				'patch14/diff_extra_2014',
				'patch14/diff_last_2014',
				'patch14/diff_plot_2014'],
		   '2.0.15' => ['patch15/diff_cast_2015',
				'patch15/diff_errout_2015',
				'patch15/diff_gnuplot_2015',
				'patch15/diff_proto_2015',
				'patch15/diff_errpari_2015',
				'patch15/diff_pari_gnuplot_2015'],
		   '2.0.16' => ['patch16/diff_gnuplot_2016'],
		   '2.1.2' =>  ['patches/diff_2.1.2_gccism'],
		  );
    if ($patches{$version}) {
      print "Patching...\n";
      foreach $patch (@{$patches{$version}}) {
	print  "cd $dir ; patch -p1 < ../$patch\n";
	system "cd $dir ; patch -p1 < ../$patch";
      }
      print "Finished patching...\n";
    }
  }

  $ftp->quit or die "Cannot quit: $!";
  $dir;
} 
