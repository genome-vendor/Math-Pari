package Math::PariBuild;

require Exporter;
@ISA = 'Exporter';
@EXPORT = qw(get_pari_version
	     pari_formatted_version
	     find_pari_dir
	     download_pari
	     patch_pari
	     download_and_patch_pari
	     make_pod
	     build_tests
	     find_paricfg
	     find_or_Configure_paricfg
	     write_paricfg
	     build_paricfg
	     find_machine_architecture
	     known_asmarch
	     inline_headers
	     is_gnu_as
	     choose_and_report_assembler
	     kernel_files
	     kernel_fill_data
	     assembler_flags
	     extra_includes
	     ep_codes_from_file
	     ep_hash_report
	     ep_in_version
	     code_C_translator
	     build_funclists
	    );

use strict;
use Config;
use File::Copy 'copy';
use File::Basename 'basename';

=head1 NAME

Math::PariBuild - utility functions used during configuration of C<Math::Pari>.

=head1 SYNOPSIS

  use Math::PariBuild;

=head1 DESCRIPTION

=over

=item C<get_pari_version($dir)>

extracts the version of GP/PARI given the build directory $dir;
version is read from $F<$dir/config/version> and has the form as in
C<"2.3.1">.  Returns undef on failure.

=cut

sub get_pari_version {
  my $dir = shift;
  my $v = "";
  open(IN, "$dir/config/version") or return;
  /(?:version|VersionMajor|VersionMinor|patch)='?(\d+(\.\d+)?)'?/
    and $v .= "$1." while <IN>;
  close(IN) or die "error closing '$dir/config/version'";
  $v =~ s/\.$// or return;
  return $v;
}

=item C<pari_formatted_version($dir)>

extracts the version of GP/PARI given the build directory $dir;
version has the form as in C<"2003001"> for version 2.3.1.  Returns
the directory name on failure.

=cut

sub pari_formatted_version {
  my $dir = shift;
  my $v;
  $v = get_pari_version $dir;
  if (defined $v) {
    $dir = $v;
  } else {
    warn(<<EOW);
Could not extract version from '$dir/config/version';
trying extract from the directory name...
EOW
  }
  return sprintf '%d%03d%03d',$1, $2, $3
    if $dir =~ /(\d+)\.(\d+).(\d+)(\.(alpha|beta))?$/;
  warn(<<EOW);
Directory `$dir' has unknown syntax...
EOW
  return $dir;
}


=item find_pari_dir()

Returns the GP/PARI build directory, looking for it as a kid, sibling,
or parent of the current directory.

=cut

my $latmus = 'src/test/in/nfields';

sub find_pari_dir {
  my ($dir, @dirs, @gooddirs);
  # Try to find alongside
  for $dir ('.', '..', '../..', '../../..') {
    @dirs = <$dir/pari-[234].*>;
    @dirs = "$dir/pari" if not @dirs and -d "$dir/pari";
    @dirs = grep -e "$_/$latmus", @dirs;
    last if @dirs;
  }
  @gooddirs = grep !/alpha|beta/, @dirs;
  @gooddirs = grep !/alpha/, @dirs unless @gooddirs;
  @gooddirs = @dirs unless @gooddirs;
  @gooddirs = sort {pari_formatted_version($a) cmp pari_formatted_version($b)}
    @gooddirs;
  return $gooddirs[-1];
}

=item download_pari()

Using FTP connection, downloads the latest version of GP/PARI, and
extracts it.  Returns the GP/PARI build directory and the version in
the format C<"2.3.1">.

=cut

sub download_pari {
  my ($srcfile) = (shift);
  my $host = 'megrez.math.u-bordeaux.fr';
  my $dir  = '/pub/pari/unix/';
  my($ftp, $ua, $base_url);

  print "Did not find GP/PARI build directory around.\n" unless defined $srcfile;

  my $match = '((?:.*\/)?pari\W*(\d+\.\d+\.\d+).*\.t(?:ar\.)?gz)$';

  my %archive;
  my $match_pari_archive = sub {
    my $file = shift;
    return unless $file =~ /$match/o;
    $file = $1;
    my $version = $2;
    if ($file =~ /alpha/) {
      $archive{alpha}{$version} = $file;
    } elsif ($file =~ /beta/) {
      $archive{beta}{$version} = $file;
    } else {
      $archive{golden}{$version} = $file;
    }
  };

  if ($srcfile and -s $srcfile) {
    die "The FILE supplied via the pari_tgz=$srcfile option did not match /$match/"
      unless $match_pari_archive->($srcfile);
  } else {
    if (-t STDIN and (-t STDOUT or -p STDOUT)) { # Interactive
      $| = 1;
      my $mess = <<EOP;
Do you want to me to fetch GP/PARI automatically?
  (If you do not, you will need to fetch it manually, and/or direct me to
   the directory with GP/PARI source via the command-line option paridir=/dir)
Make sure you have a large scrollback buffer to see the messages.
Fetch? (y/n, press Enter)
EOP
      chomp $mess;
      print "$mess ";
      my $ans = <STDIN>;
      if ($ans !~ /y/i) {
        print <<EOP;
Well, as you wish...
Rerun Makefile.PL when you fetched GP/PARI manually.
EOP
        return;
      }
    } else {
      print "Non-interactive session, autofetching...\n"
    }

    $base_url = "ftp://$host$dir";
    print "Getting GP/PARI from $base_url\n";

    eval {
      require Net::FTP;

      $ftp = Net::FTP->new($host) or die "Cannot create FTP object: $!";
      $ftp->login("anonymous","Math::Pari@")
        or die "Cannot login anonymously (",$ftp->message(),"): $!";
      $ftp->cwd($dir) or die "Cannot cwd (",$ftp->message(),"): $!";
      $ftp->binary() or die "Cannot switch to binary (",$ftp->message(),"): $!";
      my @lst = $ftp->ls();
      @lst or ($ftp->pasv() and @lst = $ftp->ls()) or die "Cannot list (",$ftp->message(),"): $!";
      #print "list = `@lst'\n";

      my $c = 0;
      %archive = ();
      for my $file (@lst) {
        $c++ if $match_pari_archive->($file);
      }
      die "Did not find any file matching /$match/ via FTP"
        unless $c;
    };
    if ($@) {
      undef $ftp;
      warn "$@\nCan't fetch file with Net::FTP, now trying with LWP::UserAgent...\n";
      # second try with LWP::UserAgent
      eval { require LWP::UserAgent; require HTML::LinkExtor }
        or die "You do not have LWP::UserAgent and/or HTML::LinkExtor installed, cannot download, exiting...";
      $ua = LWP::UserAgent->new;
      $ua->env_proxy;
      my $req = HTTP::Request->new(GET => $base_url);
      my $resp = $ua->request($req);
      $resp->is_success
        or die "Can't fetch directory listing from $base_url: " . $resp->as_string;
      my $c = 0;
      %archive = ();
      if ($resp->content_type eq 'text/html') {
        my $p = HTML::LinkExtor->new;
        $p->parse($resp->content);
        for my $link ($p->links) {
          my($tag, %attr) = @$link;
          next if $tag ne 'a';
          $c++ if $match_pari_archive->($attr{href});
        }
      } else {
        foreach my $file (split /\n/, $resp->content) {
          $c++ if $match_pari_archive->($file);
        }
      }
      die "Did not find any file matching /$match/ via FTP"
        unless $c;
    }
  }

  sub fmt_version {sprintf "%03d%03d%03d", split /\./, shift}

  my ($type, %have, %types, $best, %latest_version, %latest_file);
  for $type (qw(alpha beta golden)) {
    if ($archive{$type}) {
      $have{$type}++;
      $best = $type;
      my @files = keys %{$archive{$type}};
      print "Available $type versions: `@files'\n";
      $latest_version{$type} = (sort {fmt_version($a) cmp fmt_version($b)}
				keys %{$archive{$type}})[-1];
      $latest_file{$type} = $archive{$type}{$latest_version{$type}};
      print qq(Latest $type is `$latest_file{$type}'\n);
    }
  }

  # Special-case v2.0.14
  if (!$archive{golden} and $latest_version{beta} eq '2.0.11'
      and $latest_version{alpha} eq '2.0.14') {
    $best = 'alpha';		# It is tested!
  }

  undef $dir;
  my $version;
  if ($best) {
    my $file = $latest_file{$best};
    $version = $latest_version{$best};
    print qq(Picking $best version $version, file $file\n);
    if (-f $file) {
      print qq(Well, I already have it, using the disk copy...\n);
    } else {
      print qq(Downloading...\n);
      if ($ftp) {
        $ftp->get($file) or die "Cannot get (",$ftp->message(),"): $!";
      } else {
	my $req = HTTP::Request->new(GET => "$base_url/$file");
	my $resp = $ua->request($req);
	$resp->is_success
	  or die "Can't fetch $base_url/$file: " . $resp->as_string;
	my $base = basename($file);
	open(F, ">$base") or die "Can't write to $base: $!";
	print F $resp->content;
	close F;
      }
      print qq(Downloaded...\n);
    }
    print qq(Extracting...\n);
    my $zcat = "gzip -dc";	# zcat may be the old .Z-extractor
    print  "$zcat $file | tar -xvf -\n";
    system "$zcat $file | tar -xvf -"
      and die "Cannot extract: $!, exitcode=$?.\n";
    ($dir = $file) =~ s,(?:.*[\\/])?(.*)\.t(ar\.)?gz$,$1,
      or die "malformed name `$file'";
    -d $dir or die "Did not find directory $dir!";
  }
  if ($ftp) {
    $ftp->quit or die "Cannot quit: $!";
  }
  return ($dir, $version);
}

=item C<patches_for($version)>

Returns patches appropriate for GP/PARI version $version (formatted as in C<2.2.2>).

=cut

sub patches_for ($) {
  my ($v) = (shift);
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
		 '2.1.3' =>  ['patches/diff_2.1.3_interface'],
		 '2.1.4' =>  ['patches/diff_2.1.4_interface'],
		 '2.1.5' =>  ['patches/diff_2.1.4_interface'],
		 '2.2.2' =>  ['patches/diff_2.2.2_interface'],
		 '2.1.6' =>  ['patches/diff_2.1.6_ploth64',
			      'patches/diff_2.1.6_no-common'],
		);
  print "Looking for patches for $v...\n";
  my @p = $patches{$v} ? @{$patches{$v}} : ();
  push @p, 'patches/diff_pari-2.1.3-ix86-divl'
    if $v le '2.1.3' or $v ge '2.2' and $v le '2.2.2';
  @p;
}

=item C<patch_pari($dir [, $version])>

Applies known necessary fixes to GP/PARI build directory $dir if needed.

=cut

sub patch_pari {
  my ($dir, $version) = (shift, shift);
  $version = get_pari_version($dir) unless defined $version;
  my @patches = patches_for($version);
  if (@patches) {
    print "Patching...\n";
    my $p;
    my $patch = $Config{gnupatch} || 'patch';
    foreach $p (@patches) {
      my $cmd = "cd $dir ; $patch -p1 < ../$p";
      print "$cmd\n";
      system "$cmd"
	and warn "...Could not patch: \$?=$?, $!; continuing anyway...\n";
    }
    print "Finished patching...\n";
  }
}

=item download_and_patch_pari()

Using FTP connection, downloads the latest version of GP/PARI,
extracts it, and applies known necessary fixes if needed.  Returns the
GP/PARI build directory.

=cut

sub download_and_patch_pari {
  my ($file) = (shift);
  my ($dir, $version) = download_pari($file);
  patch_pari($dir, $version) if defined $dir;
  $dir;
}

=item C<make_pod($podfile, $gphelp_opt, $dir)>

Makes POD documentation for functions in the PARI library.  Converts
the TeX file found in GP/PARI build directory $dir to POD using the
given options for F<gphelp>.

=cut

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
  my ($targ, $how, $paridir) = @_;
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

sub scan_headers {
  my $opts = shift;
  warn "Scanning header files...\n";
  my $cmd = "$Config{cpprun} $Config{cppflags} utils/inc.h 2>&1";
  open INC, "$cmd |" or warn("Error $! from: $cmd\n"), return;
  $opts->{clk_tck_def} = 1;
  while (<INC>) {
    $opts->{have_ulong} = 1, warn "...ulong\n" if /\btypedef\b.*\bulong\s*;/;
    $opts->{clk_tck_def} = 0, warn "...CLK_TCK not defined\n"
      if /y\s*=\s*CLK_TCK\b/;
    $opts->{have_getrusage} = 1, warn "...getrusage\n" if /\bgetrusage\s*\(/;
    $opts->{have_ladd} = 1, warn "...ladd\n" if /\bladd\b/;
  }
  close INC or warn "Note (probably harmless): Errors reading from pipe: '$!', exit=$?: $cmd\n"
}

=item C<build_tests($dir)>

Converts GP/PARI test files in GP/PARI build directory $dir to Perl test suite.

=cut

sub build_tests {
  my $dir = shift;
  my $paritests = "$dir/src/test/in";
  opendir TESTS, $paritests
    or die "Cannot find tests in $paritests: $!";
  my @tests = readdir TESTS;
  closedir TESTS or die "Cannot find tests (close): $!";
  my $sou = 'test_eng/ex.t';
  my $targ = "$sou-";
  unless (-e $targ and -M $targ <= -M $sou) {
    $dir =~ s/\\/\\\\\\\\/g;
    my $quote = ($^O =~ /win32/i) ? '"' : "'";
    system "$^X -pe $quote s,CHANGE_ME,$dir, $quote $sou > $targ"
      and die "Could not run test converter: $! $?";
  }
  $sou = $targ;
  my $test;
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
}

=item C<find_paricfg($dir)>

Finds suitable (?) files F<paricfg.h> in GP/PARI build directory $dir.

=cut

sub find_paricfg {
  my $paridir = shift;
  my @paricfg = <$paridir/o.*/paricfg.h>;
  push @paricfg, <$paridir/O*/paricfg.h>;

  # Reported to work with Win32 too
  @paricfg = grep !/Odos/, @paricfg unless $^O =~ /dos|djgcc|MSWin32/i;

  # Probably not present in newer versions anymore...
  unshift @paricfg, "$paridir/win32/paricfg.h"
    if $^O eq 'MSWin32' and -f "$paridir/win32/paricfg.h";
  @paricfg;
}

=item C<find_paricfg($dir, $do_configure)>

Finds suitable (?) files F<paricfg.h> in GP/PARI build directory $dir.
If $do_configure is true, runs GP/PARI's Configure script to build
one.  Returns FALSE if F<paricfg.h> needs to be build by Perl.

=cut

sub find_or_Configure_paricfg {
  my ($paridir, $do_configure) = (shift, shift);
  my @paricfg = find_paricfg $paridir;

  return 0 unless $do_configure;
  if (@paricfg == 0) {
    print "No existing paricfg.h found, running Configure...\n";
    print  "cd $paridir ; sh ./Configure\n";
    system "cd $paridir ; sh ./Configure"
      and die "Cannot configure: $!, exitcode=$?.\n";
    print "Configuration of GP/PARI successful.\n";
    @paricfg = find_paricfg $paridir;
  }

  if (@paricfg == 0) {
    warn <<EOW;
Did not find paricfg.h.  You may need to manually copy it to libPARI
    directory from the GP/PARI build directory.
    ...Now switching to creation of paricfg.h by Perl code.
EOW
    return 0;
  }
  my $found = $paricfg[0];
  if (@paricfg > 1) {
    warn "Found multiple paricfg.h: @paricfg.\n";
    @paricfg = sort { -M $a <=> -M $b} @paricfg;
    $found = $paricfg[0];
    warn "Choosing newest paricfg.h: $found.\n";
  }
  if (-e 'libPARI/paricfg.h' and -M $found >= -M 'libPARI/paricfg.h') {
    print <<EOP;		# Duplication with build_paricfg()...
Existing libPARI/paricfg.h not older than $found.
...Will not overwrite libPARI/paricfg.h...  (remove it manually if needed);
   You may also want to remove libPARI/paricfg.h if you configuration changed
   from the time of the first build in this directory...
EOP
  } else {
    print "Found $found, copying it to libPARI...\n";
    copy $found, 'libPARI/paricfg.h'
      or die "Could not copy $found to paricfg.h: $!"
	if not -e 'libPARI/paricfg.h' or -M $found < -M 'libPARI/paricfg.h';
  }
  return 1;
}

=item write_paricfg()

Writes PARI configuration file F<libPARI/paricfg.h>.  Returns hash
with options found during the scan of the header files.

=cut

sub write_paricfg {
  my %opts;
  scan_headers(\%opts) or $opts{clk_tck_def} = 0;
  warn "Creating libPARI/paricfg.h...\n";
  open F, '> libPARI/paricfg.h' or die "open 'libPARI/paricfg.h' for write: $!";
  print F <<EOP unless $^O =~ /win32/i; # Should not we check for CygWin?
#define UNIX

EOP
  my $shellq = ($^O eq 'os2' or $^O =~ /win32/i or $^O eq 'dos') ? q(") : q(');
  print F <<EOP;
#define SHELL_Q		'\\$shellq'
EOP
  print F <<EOP;
#define GPDATADIR "/usr/local/lib/pari/galdata"
#define GPMISCDIR "/usr/local/lib/pari"

#define PARI_BYTE_ORDER    $Config{byteorder}
#define NOEXP2	/* Otherwise elliptic.t:11 rounds differetly, and fails */

EOP
  if ($opts{have_getrusage}) {
    print F <<EOP if $Config{d_times};
#define USE_GETRUSAGE 1

EOP
  } else {
    print F <<EOP if $Config{d_times} and $^O !~ /win32/i; # times() missing there...
#define USE_TIMES 1

EOP
    print F <<EOP if $Config{d_times} and $Config{i_time} and !$opts{clk_tck_def};
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

  print F <<EOP unless $opts{have_ulong};
#define ULONG_NOT_DEFINED

EOP

  my $arch = find_machine_architecture();
  my $bits64 = ($arch =~ /alpha|ia64/
		or defined($Config{longsize}) and $Config{longsize} == 8);
  print F <<EOP if $bits64;
#define LONG_IS_64BIT	1

EOP

  if (!$bits64) {
     # Order of words in a double
     my @w = unpack 'LL', pack 'd', 2;
     my $f = $w[1] ? 1 : 0;
     die "Unknown double format" unless $w[$f] == (1<<30) and $w[1-$f] == 0;
     print F <<EOP;
#define PARI_DOUBLE_FORMAT $f

EOP
  }

  print F <<EOP;
#define DL_DFLT_NAME	NULL

EOP

  print F <<EOP if $arch eq 'port';
#define __HAS_NO_ASM__

EOP

  close F or die "close 'libPARI/paricfg.h' for write: $!";
  %opts;
}

=item C<build_paricfg($dir, $do_configure)>

Builds F<libPARI/paricfg.h> either ourselves, or by looking for it in
GP/PARI build directory $dir - and running GP/PARI's Configure script
if needed.  Returns hash with options found during the scan of the
header files.

=cut

sub build_paricfg {
  my ($paridir, $do_configure) = (shift, shift);
  my %opts;
  unless (find_or_Configure_paricfg($paridir, $do_configure)) {
    # Not generated by Configure
    if (-r 'libPARI/paricfg.h') {
      print <<EOP unless $do_configure;	# Duplication with find_or_Configure_paricfg()
...Will not overwrite libPARI/paricfg.h...  (remove it manually if needed)
   You may also want to remove libPARI/paricfg.h if your configuration changed
   from the time of the first build in this directory...
EOP
    } else {
      print "...Generating libPARI/paricfg.h ...\n";
      %opts = write_paricfg();
    }
  }
  %opts;
}

# The following two functions are based on the logic in the PARI
# Configure script (updated to 2.2.8's config/arch-osname):

sub process_sparc {
  my $info = shift;
  #	    *SuperSparc*)   arch=sparcv8_super;;
  #	    *TMS390Z5[05]*) arch=sparcv8_super;; # SuperSparc I or II
  #	    *MB86934*)      arch=sparcv8_super;; # SparcLite
  #	    *RT625*)        arch=sparcv8_super;; # HyperSparc
  #	    *CY605*)        arch=sparcv8_super;;
  return 'sparcv8_super' if $info =~ /SuperSparc|TMS390Z5[05]|CY605|MB86934|RT625/;

  #	    *TMS390S1[05]*) arch=sparcv8_micro;; # MicroSparc I
  #	    *MB86904*)      arch=sparcv8_micro;; # MicroSparc II
  #	    *MB86907*)      arch=sparcv8_micro;; # TurboSparc
  return 'sparcv8_micro' if $info =~ /TMS390S1[05]|MB8690[47]/;
  return shift;
}

=item find_machine_architecture()

Returns the type of the processor of the current machine.

=cut

sub find_machine_architecture () {
  my $os = (split ' ', $Config{myuname})[0];

  my $machine = $os;		# Handles fx2800
  if ($os =~ /^irix/) {
    $machine = 'irix';
  } elsif ($os =~ /^hp/) {
    $machine = `uname -m` || 'hppa';
  } elsif ($os eq 'os2' or $os eq 'netbsd'
	   or $os eq 'freebsd' or $os =~ /^cygwin/) {
    chomp($machine = `uname -m`);
    $machine ||= 'ix86';
  } elsif (0 and $os =~ /win32/i and not $Config{gccversion}) {
    # Not needed with rename of kernel1.s to kernel1.c?
    $machine = 'port'; # Win32 compilers would not understand the assmebler anyway
  } elsif ($os eq 'ultrix') {
    $machine = 'mips';
  } elsif ($os eq 'nextstep' or -d '/NextApps') {
    chomp($machine = `file /bin/sh | sed 's/.*(for architecture \(.*\))/\1/'`);
  } elsif ($os eq 'osf1') {
    $machine = 'alpha' if (split ' ', $Config{myuname})[4] eq 'alpha';
  } elsif ($os =~ /^cygwin/) {
    $machine = $ENV{HOSTTYPE};
  } elsif ($os eq 'linux') {
    chomp($machine = `uname -m`);
    $machine = 'sparcv9' if $machine eq 'sparc64';
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
    } elsif ($type =~ /^sun4[ce]?/) {
      $machine = 'sparcv7';
    } elsif ($type =~ /^sun4[dm]/) {
      local $ENV{PATH} = "$ENV{PATH}:/dev/sbin";
      my $info = `(prtconf||devinfo)2>&-`;
      $info = join ' ', grep /(TI|FMI|Cypress|Ross),/, split "\n", $info;
      $machine = process_sparc $info, 'sparcv8';
    } elsif ($type eq 'sun4u') {
      $machine = 'sparcv9';
    } elsif ($type =~ /^i.*pc$/) {
      $machine = 'ix86';
    } elsif ((split ' ', $Config{myuname})[3] eq 'sun') {
      $machine = 'm86k';
    }
  } elsif ($os eq 'gnu') {
    chomp($machine = `uname -m`);
    $machine = 'ix86' if $machine =~ /^i386-/;
  }

  if ( $machine ne 'alpha'
       and defined($Config{longsize}) and $Config{longsize} == 8 ) {
    $machine = 'port';			# No assembler for 64bit - unless alpha
  }

  # For older PARI:
  ###  $machine = 'sparcv8super'
  ###    if $machine eq 'sparcv9' or $machine eq 'sparcv8_hyper'
  ###       or $machine eq 'sparcv8_super';
  ###  $machine = 'sparcv8micro' if $machine eq 'sparcv8_micro';

  # This part is probably not needed and never entered
  if (not defined $machine
      and $Config{myuname}
      =~ /\b(sun3|sparcv7|sparcv8_micro|sparcv8_super|alpha|hppa|[ix]\d86)\b/) {
    $machine = $1;
  } elsif (not defined $machine) {
    chomp($machine = `uname -m`);
  }
  $machine =~ s/[ix]\d86/ix86/ if defined $machine;

  print "...Processor of family `$machine' detected\n";
  return $machine;
}

# Which files to catenate to produce pariinl.h.  Apparently, the only
# need to go to asm0.h is to undo the effect of ASMINLINE in paricfg.h.
# Note that we do ASMINLINE from the command line.

sub sparcv8_inl {
  my ($asmarch, $pari_version) = (shift, shift);
  return ['none/asm0.h','none/level1.h']
    if $Config{osname} =~ /^(linux|nextstep)$/;
  return ['sparcv8/level0.h','none/level1.h'] if $pari_version < 2002006;
  return ['sparcv8_micro/level0_common.h','sparcv8_micro/level0.h',
	  'none/level1.h'] if $asmarch eq 'sparcv8_micro';
  return ['sparcv8_micro/level0_common.h','none/divll.h',
	  'none/level1.h'] if $asmarch eq 'sparcv8_super';
  # No for sparcv8...
}

sub inline_headers_arr {     # These files are cat()ed to pariinl.h
  my ($asmarch, $pari_version) = (shift, shift);
  return sparcv8_inl($asmarch, $pari_version) if $asmarch =~ /^sparcv8/;
  my %h = (
	       alpha	      => ['none/asm0.h','none/level1.h'],
	       hppa	      => ['none/asm0.h','none/level1.h'],
	       ix86	      => ['ix86/level0.h','none/level1.h'],
	       m86k	      => ['none/level0.h','none/level1.h'],
	       none	      => ['none/level0.h','none/level1.h'],
	       # ppc is not done yet (2.0.15)
	($pari_version > 2002007
		? (ppc		     => ['ppc/asm0.h', 'none/divll.h'],
		   ia64		     => ['ia64/asm0.h','ia64/asm1.h'])
		: ()),
	       sparcv7	      => ['none/asm0.h','none/level1.h'],
#	       sparcv8	      => $sparcv8_inl,
#	       sparcv8_micro  => $sparcv8_inl,
#	       sparcv8_super  => $sparcv8_inl,
	       # sparcv9 is not done yet (2.0.15)
  );
  $h{$asmarch};
}

sub inline_headers {
  my ($asmarch, $pari_version) = (shift, shift);
  my $inlines = inline_headers_arr($asmarch, $pari_version)
	or die "Unknown inlines for '$asmarch'";
  my @inlines = @$inlines;
  unshift @inlines, 'none/int.h' if $pari_version >= 2002005;
  unshift @inlines, 'none/tune.h' if $pari_version >= 2002008;
  map "\$(PARI_DIR)/src/kernel/$_", @inlines;
}

sub known_asmarch {
  defined inline_headers_arr(@_);
}

sub is_gnu_as {
  local $/;
  my $ass = $ENV{AS} || 'as';
  open ASS, "$ass --version 2>&1 |";
  my $assout;
  eval {
    local $SIG{ALRM} = sub {die};
    eval {alarm 10};			# Be extra safe...
    $assout = <ASS>;
    close ASS;
    unless ($assout) {
      open ASS, "$ass -v 2>&1 |";
      $assout = <ASS>;
      close ASS;
    }
    eval {alarm 0};
  };
  return ($assout and $assout =~ /GNU/);
}

sub choose_and_report_assembler {
  my $machine = shift;
  my %asmarch = (
	      sun3	   => 'm86k',
	      sparc	   => 'sparcv8_micro',
	      sparcv9	   => 'sparcv8_micro',
	      port	   => 'none',
	      mips	   => 'none',
	      fx2800	   => 'none',
	      hppa	   => ($Config{osvers} =~ /^.\.10\./
			       ? 'hppa' : 'none'),
	     );
  my $asmarch = $asmarch{$machine} || $machine; # Temporary only
  unless (known_asmarch $asmarch) {
    warn <<EOW;
####	Do not know how to build for assembler `$asmarch'.	####
####	Reversing to assembler-less type `port'.		####
####						 		####
####	If you think your processor's assembler is supported	####
####	by PARI, edit libPARI/Makefile.PL and report.		####
####						 		####
####	Alternatively, specify machine=YOURTYPE on the		####
####	  perl Makefile.PL line			 		####
####	Recognized types:			 		####
####	  alpha hppa m86k none sparcv7 sparcv8 sparcv8_micro	####
####	  sparcv8_super ix86			 		####
EOW
    $machine = 'port';
    $asmarch = 'none';
  }
  if ($asmarch eq 'none') {
    print "...I will  use portable assembler-less build\n";
  } else {
    print "...I will use assembler build of type '$asmarch'.\n";
  }

  print <<EOP if $asmarch eq 'hppa';
###
###  Some time ago HPPA assembler files were not relocatable,
###    if this is still true, they are probably unsuitable for dynamic linking.
###  It is advisable to restart Makefile.PL with an extra argument
###     machine=port
###  if you are planning for dynamic linking of Math::Pari.
###
###  NOTE:  machine=port results in a significant drop in performance.
###  For a static build (which makes a new perl executable with the library
###	  compiled in [and arranges for it to be compiled in when
###	  other extensions are statically built later]):
###    perl Makefile.PL LINKTYPE=static
###    make static
###    make perl
###    make test
###    make install
###
EOP
  return $asmarch;
}

# Output (last two optional):
#   [Which file to compile, whether you need to preprocess it to ./kernel1.s,
#    Additional file to compile, need? to preprocess it to ./kernel2.s,]

sub sparcv8_kernel_files_old {
  my ($asmarch, $pari_version, $Using_gnu_as) = (shift, shift, shift);
  my $_ext = (($pari_version < 2000015) ? 's' : 'S');
  my $sparcv8_cvt = $Using_gnu_as || $Config{osname} =~ /^(linux|nextstep)$/;
  my $sparcv8_kernel = ($sparcv8_cvt
			? ["sparcv8/level0.$_ext", 1,
			   "sparcv8/level0_$asmarch.$_ext", 1]
			: ["sparcv8/level0.$_ext", 0,
			   "sparcv8/level0_$asmarch.$_ext", 0]);

  # kernel2.o is not needed if the compiler can inline assembler:
  my $sparcv8_need_kernel2 =
    !$Config{gccversion} || $Config{osname} =~ /^(linux|nextstep)$/;
  $sparcv8_kernel = [$sparcv8_kernel->[2], $sparcv8_kernel->[3]]
    unless $sparcv8_need_kernel2;
  return   $sparcv8_kernel;
}

sub kernel_files {
  my ($asmarch, $pari_version, $Using_gnu_as) = (shift, shift, shift);

  return sparcv8_kernel_files_old($asmarch, $pari_version, $Using_gnu_as)
    if $asmarch =~ /^sparcv8/ and $pari_version < 2002006;

  my $sparcv8_kernel = ["sparcv8_micro/level0_common.S", 1,
			"$asmarch/level0.S", 1];	# 2.2.* only

  # Default ["$asmarch/level0.s", 0]
  my %level0 = (
		alpha		     => '',
		hppa		     => '', # was ['none/level0.c', 0] before 2.0015
	($pari_version > 2002007
		? (ppc		     => ["none/level0.c", 0],
		   ia64		     => ["ia64/level0.s", 0])
		: ()),
		ix86		     => ['ix86/l0asm.c',  1],
		m86k		     => ["none/level0.c", 0],
		none		     => ["none/level0.c", 0],
		# ppc is not done yet (2.0.15)
		sun3		     => '',
		sparcv7		     => '',
#		sparcv8		     => $sparcv8_kernel,
		sparcv8_micro	     => $sparcv8_kernel,
		sparcv8_super	     => $sparcv8_kernel,
		# sparcv9 is not done yet (2.0.15)
	       );

  return $level0{$asmarch} || ["$asmarch/level0.s", 0];
}

sub kernel_fill_data {
  my ($kernels, $hash) = (shift, shift);
  # The original file
  $hash->{file1}	= "\$(PARI_DIR)/src/kernel/$kernels->[0]";
  $hash->{convert1}	= $kernels->[1];

  # The (possibly) preprocessed file
  $hash->{converted1}	= $hash->{convert1} ? 'kernel1.s' : $hash->{file1};

  # Extra bookkeeping
  ($hash->{header1}	= $hash->{file1}) =~ s/\.c$/.h/;
  $hash->{header1}	= '' unless -r $hash->{header1}; # Put as a dependence
  ($hash->{dir1}	= $hash->{file1}) =~ s/[^\/]*\.[csS]$//; # Include with -I

  # Additional file to compile (original and converted)
  $hash->{file2}	= $hash->{converted2} = "\$(PARI_DIR)/src/kernel/$kernels->[2]";
  $hash->{file2}	= $hash->{converted2} = '' unless $kernels->[2];
  $hash->{converted2}	= 'kernel2.s' if $kernels->[3];

  $hash->{convert}	=  ($hash->{converted2} eq 'kernel2.s'
			    or $hash->{converted1} eq 'kernel1.s');
  if ( $^O =~ /win32/i
       and $Config{cc} =~ /\bcl/ # M$ VC doesn't understand .s
       and $hash->{converted1} eq 'kernel1.s' ) {
    $hash->{converted1} = 'kernel1.c';
  }

}

sub assembler_flags {
  my ($machine, $Using_gnu_as) = (shift, shift);
  my %assf  = (
	       # alpha  => "-O1",		# Not supported any more
	       sparc  => ($Config{osname} eq 'solaris'
			  ? "-P -T -I."
			  : "-P -I."),
	       hppa   => "+DA1.1",
	      );
  my $assflags  = $assf{$machine =~ /sun3|sparc/ ? 'sparc' : $machine} || '';
  $assflags .= ' -D__GNUC__'	# hiremainder problem with gcc on Solaris
    if not $Using_gnu_as and $Config{gccversion};
  return $assflags;
}

sub extra_includes {
  my $pari_dir = shift;
  # Some #include directives assume us inside $pari_dir/OARCH; replace by src
  return join ' -I ', '', grep -d, "$pari_dir/src/systems/$^O", "$pari_dir/src";
}

sub build_funclists {
  my $pari_dir = shift;
  return unless -d "$pari_dir/src/desc"; # Old version, no autogeneration
  return if -f "$pari_dir/src/language/init.h"
        and -f "$pari_dir/src/desc/pari.desc";
  open FL, "> $pari_dir/src/funclist" and close FL	# Ignore errors
    unless -f "$pari_dir/src/funclist";
  (system("cd $pari_dir/src/desc && make")
     and system("cd $pari_dir/src/desc && make SHELL=cmd")
   or not -s "$pari_dir/src/desc/pari.desc") and
      (unlink("$pari_dir/src/desc/pari.desc"),
       die <<EOW);
###
###  Apparently, we failed to build function descriptions of GP/PARI.
###  Try editing $pari_dir/src/desc/Makefile - a typical reason
###  is a wrong value of SHELL for your system.  You can run make in
###  $pari_dir/src/desc manually too...
EOW
}

=item ep_codes_from_file($filename,%hash,%names)

Adds to the %hash the string interface descriptions corresponding to
the numeric codes use in the file's entree array.  %hash is indexed by
by the numeric codes; the value are references to arrays with the corresponding
string interface descriptions.

Adds to %names the list of name => code values.

=cut

sub ep_codes_from_file ($\%\%) {
  my ($file, $descrh, $names) = (shift, shift, shift);
  local $_;
  open IN, "< $file" or warn "Cannot open `$file': $!" and return;
  while (<IN>) {
    next unless /^\s*{\s*\"/;
    chomp;
    warn("Unrecognized line: `$_'\n"), next
      unless /^\s*\{\s*"(\w+)"\s*,\s*(\d+)\s*,[^,]*,\s*\d+\s*,(?:\s*\d+\s*,)?\s*("((?:\\.|[^"])*)"|NULL)\s*(,|\})/;
    next unless defined $4;
    #print;
    my ($name, $code, $descr) = ($1, $2, $4);
    $descrh->{$code} = [] unless exists $descrh->{$code};
    push @{$descrh->{$code}}, $descr
      unless grep $descr eq $_, @{$descrh->{$code}};
    warn "! Duplicate code $code for function '$name' (was $names->{$name})\n"
      if defined $names->{$name};
    $names->{$name} = [$code, $descr];
  }
  close IN or warn "Cannot close `$file': $!";
}

=item ep_hash_report(%hash, %names, $fh)

Writes to $fh the diagnostic about problemes with the string interface
descriptions corresponding to the numeric codes.  If $fh is false,
returns TRUE if no problem were found.

=cut


my $expected_codes_as_in = '2.1.3';
# Not in 2.1.3
my ($dummy1, %old_expected_codes) = split /\s+/, qq(
			13	GD0,L,D0,G,
			34	vLLL
	);
my ($dummy2, %expected_codes) = split /\s+/, qq(
			1	Gp
			2	GG
			3	GGG
			4	GGGG
			10	lG
			11	L
			12	GnP
			14	GDn
			16	ls
			18	G
			19	vLL
			20	lGG
			21	GL
			22	GVI
			23	GL
			24	LG
			25	GGD0,L,
			26	GnG
			27	V=GIp
			28	GDVDI
			29	GGp
			30	lGGG
			31	GDGDGD&
			32	GGL
			33	GGGD0,L,p
			35	vLGG
			37	V=GGIp
			45	LGD0,L,
			47	V=GGIDG
			48	V=GGIDG
			49	GGDVDVDI
			57	vLs
			59	vLGGGG
			62	GD0,G,D0,G,D0,L,p
			73	LV=GGIpD0,L,D0,L,
			83	vV=GGI
			84	vGVI
			85	vS
			86	vV=GGGI
			87	vV=GID0,L,
			91	GD0,L,DGp
			96	GD0,L,DGp
		       );
# Some historic changes in interfaces we do not care about (E vs I)
my $t;
my %variations = map {($t = $expected_codes{$_}) =~ s/I/E/ ? ($_, $t) : ()}
  keys %expected_codes;

my %known_unimplemented = (57 => 1, 62 => 1);

sub ep_hash_report (\%;\%$) {
  my ($h, $names, $fh) = (shift, shift, shift);
  my ($c, @list);
  $names = {} unless defined $names;
  my @keys = grep {$_ ne 0 and $_ ne 99} keys %$h;
  if (@list = grep {not exists $h->{$_}} keys %expected_codes) {
    return unless $fh;
    print $fh <<EOP;
  Cosmetic only: the following numeric interfaces are not used any more:
EOP
    for $c (sort @list) {
      print $fh <<EOP;
    $c (was meaning "$expected_codes{$c}" in $expected_codes_as_in)
EOP
    }
  }
  if (@list = grep {not exists $expected_codes{$_}} @keys) {
    return unless $fh;
    print $fh <<EOP;
  Harmless: the following numeric interfaces are new:
EOP
    for $c (sort @list) {
      print $fh <<EOP;
    $c meaning "@{$h->{$c}}"
EOP
      my $list = join ", ", grep $names->{$_}[0] == $c, sort keys %$names;
      print $fh <<EOP;
      (appears for $list)
EOP
    }
  }
  if (@list = grep @{$h->{$_}} != 1, @keys) {
    return unless $fh;
    print $fh <<EOP;
  May be harmless: non-unique string interfaces for numeric interfaces:
EOP
    for $c (sort @list) {
      print $fh <<EOP;
    $c meaning "@{$h->{$c}}"
EOP
      my $list = join ", ", grep $names->{$_}[0] == $c, sort keys %$names;
      print $fh <<EOP;
      (appears for $list)
EOP
    }
  }
  if (@list = grep {exists $expected_codes{$_}
		      and "@{$h->{$_}}" ne $expected_codes{$_}
			and "@{$h->{$_}}" ne $variations{$_}
			  and not $known_unimplemented{$_}} @keys) {
    return unless $fh;
    print $fh <<EOP;
  Possible problems with Math::Pari interface to GP/PARI:
	unexpected value of string interfaces for numeric interfaces:
EOP
    for $c (sort @list) {
      print $fh <<EOP;
    $c meaning "@{$h->{$c}}" (was meaning "$expected_codes{$c}" in $expected_codes_as_in)
EOP
      my $list = join ", ",
	grep {$names->{$_}[0] == $c and $names->{$_}[1] ne $expected_codes{$c}}
	  sort keys %$names;
      print $fh <<EOP;
      (may affect functions $list)
EOP
    }
  }
  return 1;
}

=item ep_in_version($version)

Updates the list of codes for the given version of GP/PARI (formatted as in
2002002).

=cut

sub ep_in_version ($) {
  my $v = shift;
  if ($v >= 2002002) {
    my $c;

    for $c (qw(26 62)) {
      delete $expected_codes{$c};
    }
  }
}

=item code_C_translator()

Returns string for C code to translate code string to the interface number.

Due to a bug in C_constant(), need to translate C<''> to 9900 by hand outside
of this subroutine.

=cut

sub code_C_translator {
  # Some historic changes in interfaces we do not care about (E vs I)
  my %c = (%old_expected_codes, %expected_codes);
  my %codes;
  @codes{values %c} = keys %c;
  my $k;
  for $k (keys %codes) {
    (my $kk = $k) =~ s/I/E/g;
    $codes{$kk} = $codes{$k} unless exists $codes{$kk};
    ($kk = $k) =~ s/D0,G,/DG/g;			# New alternative syntax
    $codes{$kk} = $codes{$k} unless exists $codes{$kk};
  }
  #$codes{''} = 9900;		# bug in C_constant - can't handle
  $codes{'p'} = 0;

  require ExtUtils::Constant;
  my @t = ExtUtils::Constant::constant_types(); # macro defs
  my @tt =
    ExtUtils::Constant::C_constant(
	  'Math::Pari::func_type', 
          'func_ord_by_type', undef, undef, undef, undef,
          map {{name => $_, value => $codes{$_}, macro => 1}} keys %codes);

  # 23 == 21  47==48  91 == 96 (this one unsupported)
  join '', @t, @tt;
}

=back

=cut

1;

