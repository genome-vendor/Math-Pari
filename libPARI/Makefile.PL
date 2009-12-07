use ExtUtils::MakeMaker;
use Config '%Config';
# Based on Makemakefile shell script by Bruno Haible 10.12.1992,
# modified by Babecool 23.7.1993
# Created by Ilya Zakharevich <ilya@math.ohio-state.edu>

$main_paridir = $common::main_paridir;
die "Please start Makefile.PL at toplevel, not in libPARI"
  unless $main_paridir;

if ($main_paridir =~ m,^([a-z]:)?/,i) {
  $main_paridir = $common::main_paridir;
} else {
  $main_paridir = "../$common::main_paridir";
}

$pari_version = $common::pari_version;

# The following two functions are based on the logic in the PARI
# Configure script:

sub process_sparc {
  my $info = shift;
#	    *TMS390Z5[05]*) arch=sparcv8_super;; # SuperSparc I or II
#	    *CY605*)        arch=sparcv8_super;;
#	    *MB86934*)      arch=sparcv8_super;; # SparcLite
#	    *RT625*)        arch=sparcv8_super;; # HyperSparc
  return 'sparcv8_super' if $info =~ /TMS390Z5[05]|CY605|MB86934|RT625/;
#	    *TMS390S1[05]*) arch=sparcv8_micro;; # MicroSparc I
#	    *MB86904*)      arch=sparcv8_micro;; # MicroSparc II
#	    *MB86907*)      arch=sparcv8_micro;; # TurboSparc
  return 'sparcv8_micro' if $info =~ /TMS390S1[05]|MB8690[47]/;
  return shift;
}

sub process_machine {
  my $os = (split ' ', $Config{myuname})[0];
  my $machine = $os;
  if ($os =~ /^hp/) {
    $machine = 'hppa';
  } elsif ($os eq 'os2' or $os eq 'freebsd') {
    $machine = 'ix86';
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
      $machine = process_sparc $info, 'sparcv8';
    } elsif ($type eq 'sun4u') {
      $machine = 'sparcv9';
    } elsif ($type =~ /^i.*pc$/) {
      $machine = 'ix86';
    } elsif ((split ' ', $Config{myuname})[3] eq 'sun') {
      $machine = 'm86k';
    }
  }
  # For older PARI:
###  $machine = 'sparcv8super'
###    if $machine eq 'sparcv9' or $machine eq 'sparcv8_hyper'
###       or $machine eq 'sparcv8_super';
###  $machine = 'sparcv8micro' if $machine eq 'sparcv8_micro';

  # This part is probably not needed and never entered
  if (not defined $machine
      and $Config{myuname}
      =~ /\b(sun3|sparcv7|sparcv8_micro|sparcv8_super|alpha|hppa|[ix]\d86)\b/)
  {
    $machine = $1;
  } elsif (not defined $machine) {
    chomp($machine = `uname -m`);
  }
  $machine =~ s/[ix]\d86/ix86/ if defined $machine;

  print "...Processor of family `$machine' detected\n";
  return $machine;
}

# This part (same as the functions above) follows Configure of PARI.

if (not defined $machine and defined $common::machine) {
  $machine = $common::machine;
} elsif (not defined $machine and @ARGV and $ARGV[0] =~ /^machine=(.*)$/i) {
  $machine = $1; shift;
}

if (defined $machine) {
  print "...Via command-line: processor family `$machine'\n";
} else {
  $machine = process_machine;
}

%asmarch = (
	     sun3	   => 'm86k',
	     sparc	   => 'sparcv8_micro',
	     sparcv9	   => 'sparcv8_micro',
	     port	   => 'none',
	     mips	   => 'none',
	     fx2800	   => 'none',
	     hppa	   => ($Config{osvers} =~ /^.\.09/
			       ? 'none' : 'hppa'),
	    );

$asmarch = $asmarch{$machine} || $machine;
if ($asmarch eq 'none') {
  print "...I will use portable assembler-less build\n";
} else {
  print "...I will use assembler build of type '$asmarch'.\n";
}

print <<EOP if $asmarch eq 'hppa';
###
###  Apparently current HPPA assembler files are not relocatable,
###  thus probably unsuitable for dynamic linking.
###  It is advisable to restart Makefile.PL with an extra argument
###     machine=port
###  if you are planning for dynamic linking of Math::Pari.
###
###  NOTE:  machine=port results in a significant drop in performance.
###  For a static build:
###    perl Makefile.PL LINKTYPE=static
###    make static
###    make perl
###    make test
###    make install
###
EOP

# This part is based on analysing Makefile.SH in */kernel/* subdirectories.

$sparcv8_have_kernel2 = 
  !$Config{gccversion} || $Config{osname} =~ /^(linux|nextstep)$/;

# Using gnu assembler XXXX: how to find?

if ($machine =~ /sparc|alpha|hppa/) {
  $asscmd = '$(AS) $(ASFLAGS)';	#  $(CCCDLFLAGS)?
} else {
  $asscmd = '$(CCCMD) $(CCCDLFLAGS)';
}

{
  local $/;
  $ass = $ENV{AS} || 'as';
  open ASS, "$ass --version 2>&1 |";
  $assout = <ASS>;
  close ASS;
  unless ($assout) {
    open ASS, "$ass -v 2>&1 |";
    $assout = <ASS>;
    close ASS;
  }
  $Using_gnu_as = 1 if $assout and $assout =~ /GNU/;
}

$sparcv8_cvt = $Using_gnu_as || $Config{osname} =~ /^(linux|nextstep)$/;
$sparcv8_inl = ( $sparcv8_cvt
	       ? ['none/asm0.h','none/level1.h']
	       : ['sparcv8/level0.h','none/level1.h'] );

# These files are cat()ed to pariinl.h
%inlines = (
	     alpha	   => ['none/asm0.h','none/level1.h'],
	     hppa	   => ['none/asm0.h','none/level1.h'],
	     ix86	   => ['ix86/level0.h','none/level1.h'],
	     m86k	   => ['ix86/level0.h','none/level1.h'],
	     none	   => ['none/level0.h','none/level1.h'],
	    # ppc is not done yet (2.0.15)
	     sparcv7	   => ['none/asm0.h','none/level1.h'],
	     sparcv8	   => $sparcv8_inl,
	     sparcv8_micro => $sparcv8_inl,
	     sparcv8_super => $sparcv8_inl,
	    # sparcv9 is not done yet (2.0.15)
	   );

unless (exists $inlines{$asmarch}) {
  warn <<EOW;
####	Do not know how to build for assembler `$asmarch'.	####
####	Reversing to assembler-less type `port'.		####
####						 		####
####	If you think your processor's assembler is supported	####
####	by PARI, edit libPARI/Makefile.PL and report.		####
EOW
  $machine = 'port';
  $asmarch = 'none';
}

$inlines = $inlines{$asmarch} or die "Unknown inlines for '$asmarch'";
@inlines = map "$main_paridir/src/kernel/$_", @$inlines;

# Based on what to do with kernel$(OBJ_EXT).

# Format (last two optional):
#   [Which file to compile, whether you need to preprocess it to ./kernel1.s,
#    Additional file to compile, need? to preprocess it to ./kernel2.s,]

$_ext = (($pari_version < 2000015) ? 's' : 'S');

$sparcv8_kernel = ($sparcv8_cvt
		  ? ["sparcv8/level0.$_ext", 1,
		     "sparcv8/level0_$asmarch.$_ext", 1]
		  : ["sparcv8/level0.$_ext", 0,
		     "sparcv8/level0_$asmarch.$_ext", 0]);

$sparcv8_kernel = [$sparcv8_kernel->[0], $sparcv8_kernel->[1]]
  unless $sparcv8_have_kernel2;

# Default ["$asmarch/level0.s", 0]
%level0 = (
	     alpha	   => '',
	     hppa	   => '', # was ['none/level0.c', 0] before 2.001500
	     ix86	   => ['ix86/l0asm.c', 1],
	     m86k	   => ["none/level0.c", 0],
	     none	   => ["none/level0.c", 0],
	    # ppc is not done yet (2.0.15)
	     sun3	   => '',
	     sparcv7	   => '',
	   # Following 2 guys should be fixed on linux/NS
	     sparcv8	   => $sparcv8_kernel,
	     sparcv8_micro => $sparcv8_kernel,
	     sparcv8_super => $sparcv8_kernel,
	    # sparcv9 is not done yet (2.0.15)
	    );

$kernel1machine = $level0{$asmarch} || ["$asmarch/level0.s", 0];

# The original file
$kernel1		= "$main_paridir/src/kernel/$kernel1machine->[0]";
$kernel1convert	= $kernel1machine->[1];

# The compiled file
$kernel1c	= $kernel1convert ? 'kernel1.s' : $kernel1;

# Extra bookkeeping
($kernel1h	= $kernel1) =~ s/\.c$/.h/;
$kernel1h	= '' unless -r $kernel1h;
($kernel1dir	= $kernel1) =~ s/[^\/]*\.c$//;

# Additional file to compile (original and converted)
$kernel2	= $kernel2c = "$main_paridir/src/kernel/$kernel1machine->[2]";
$kernel2	= $kernel2c = '' unless $kernel1machine->[2];
$kernel2c	= 'kernel2.s'
  if $kernel1machine->[3];

%mp = ( );
$mp_machine = $mp{$asmarch} || 'none/mp.c';
$mp = "$main_paridir/src/kernel/$mp_machine";

%mpinl = ( );
$mpinl_machine = $mpinl{$asmarch} || 'none/level1.c';
$mpinl = "$main_paridir/src/kernel/$mpinl_machine";
($mpinlh = $mpinl) =~ s/\.c$/.h/;

$libs="-lm";
$sol_under = ''; 

#@cfiles = ();
#opendir DIR, '.' or die "Cannot opendir . : $!";
#while (defined ($_ = readdir DIR)) {
#  push @cfiles, $_ if /\.c$/i 
#    and not /^(plot|ix86|version|mpin|dummy)/;
#}
#closedir DIR or die "Cannot closedir . : $!";

@cfiles = <$main_paridir/src/*/*.c>;
@cfiles = grep !/\/(plot|ix86|version|mpin|dummy|gp_rl|gp\.)/, @cfiles;
%cfiles = map {m,/([^/.]*)\.c,i; ($1,$_)} @cfiles;


if ($Config{osname} eq 'solaris') {
  @sc_dirs = </opt/SUNWspro/SC*>;
  die "Cannot find SUNWspro dirs, needed for -lsunmath" unless @sc_dirs;
  $libs .= " -L$sc_dirs[-1]/lib -lsunmath";
  $sol_under = '_nounder';	# Need to filter
} elsif ($Config{osname} eq 'sunos') {
  # $machine = 'sparcv8_super' if $Config{myuname} =~ /\bsun4m\b/o;
  $mycflags .= ' -DULONG_NOT_DEFINED';
}

# We remove optimize options, since Perl probably knows them
# and one can specify "OPTIMIZE=-O2 -m486" on the command line

if ($Config{cc} eq 'gcc') {
  #$mycflags="-O2 -g";		# We do not put -ansi to avoid versioncflags
  $mycflags="";		# We do not put -ansi to avoid versioncflags
} elsif ($machine eq 'hppa') {
  #$mycflags="-O -Aa -DHPPA";
  $mycflags="-Aa -DHPPA";
} elsif ($machine eq 'ix86') {
  #$mycflags="-O2 -m486";
  $mycflags="";
} else {
  #$mycflags="-O";
  $mycflags="";
}

# These are not needed with newer versions of PARI
# $mycflags .= ' -Dshifts=pari_shifts'; # Conflict with libnsl

$my_redefines = '-Derr=pari_err';
$mycflags .= " $my_redefines"; # On linux it can get a wrong dynamic loading

$noexp2 = '';

if ($machine eq 'alpha') {
  $mycflags .= " -DLONG_IS_64BIT";
} elsif ($Config{osname} eq 'solaris') {
  $mycflags .= " -DSOLARIS";
} elsif ($Config{osname} eq 'os2') {
  $noexp2 = 1;
  $mycflags .= ' -DMALLOC_PROCS -DHAS_STRICMP';
  $add_ar_flags = '-p64';	# Need big library. 32 is OK without debugging
} elsif ($Config{osname} eq 'linux') {
  $noexp2 = 1;
}

$mp_cflags = '';
if ($machine eq 'sun3') {
  $mp_cflags = '\$(CPlPl) \$(NOASM)'; # NOASM seems to be empty
} elsif ($machine eq 'port') {
  push @portfiles, 'mp.c';
  $mycflags .= " -D__HAS_NO_ASM__";
}

$cfiles{plotgnuplot} = "$main_paridir/src/graph/plotgnuplot.c";
$cfiles{plotport} = "$main_paridir/src/graph/plotport.c";

@obj_files = map { "$_\$(OBJ_EXT)" } keys %cfiles;
push @obj_files, 'kernel$(OBJ_EXT)', 'mp$(OBJ_EXT)', 'mpinl$(OBJ_EXT)',
  ($kernel2c
   ? 'kernel2$(OBJ_EXT)'
   : ());

#unless (-f 'genpari.h') {
#  $genpari = ($machine eq 'sun3' ? 'genpari68k.h' : 'genpariother.h');
#  eval {link $genpari, 'genpari.h'} or do {
#    require File::Copy; 
#    File::Copy::copy ($genpari, 'genpari.h') 
#      or die "Cannot copy $genpari to genpari.h";
#  }
#}

$mycflags .= " -DDYNAMIC_PLOTTING";
# OMF build needs no underscores:
$mycflags .= " -D__NO_AOUT" if $^O eq 'os2' and !$OS2::is_aout;

%assflags = (
	     alpha  => "-O1",
	     sparc  => ($Config{osname} eq 'solaris' 
			? "-P -T -I."
			: "-P -I."),
	     hppa   => "+DA1.1",
	    );

$assflags{$machine} = $assflags{sparc} if $machine =~ /sun3|sparc/;

WriteMakefile(
    NAME	=> 'Math::PARI::libPARI',
    LINKTYPE	=> 'static',
    LIBS	=> $libs,
    OBJECT	=> join(' ', @obj_files), 
    macro	=> { 
		    MP_CFLAGS	  => $mp_cflags,
		    SOL_UNDER	  => $sol_under, 
		    ADD_AR_OPT	  => ($add_ar_flags || ''),
		    CPP		  => $Config{cppstdin},
		    CPPMINUS	  => $Config{cppminus},
		    ASSCMD	  => "$asscmd $my_redefines",
		    ASFLAGS	  => $assflags{$machine} || '',
		   },
    DEFINE	=> $mycflags,
    INC		=> "-I$main_paridir/src/headers -I$main_paridir/src/graph -I.",
    C		=> \@cfiles,
    SKIP	=> [qw( distclean test dist makeaperl xs_o static)],
    clean	=> {'FILES' => 
		    'libPARI$(LIB_EXT) pariinl.h kernel1.s kernel2.s'},
);

sub MY::top_targets {
  my $converted = '';
  $converted = <<EOC if $kernel2c eq 'kernel2.s' or $kernel1c eq 'kernel1.s';
kernel1.s: $kernel1 $kernel1h
	\$(CPP) -I. \$(INC) -I$kernel1dir \$(DEFINE) \$(CPPMINUS) < $kernel1 | perl -ne "s/%\\s+/%/g; print unless /^\\s*#/" > \$@

kernel2.s: $kernel2
	\$(CPP) -I. \$(INC) -I$kernel1dir \$(DEFINE) \$(CPPMINUS) < $kernel2 | perl -ne "s/%\\s+/%/g; print unless /^\\s*#/" > \$@

EOC
  '
all :: libPARI$(LIB_EXT)

static ::       libPARI$(LIB_EXT)

test:

libPARI$(LIB_EXT): $(OBJECT)
	-$(RM_F) libPARI$(LIB_EXT)
	$(AR) $(ADD_AR_OPT) cr libPARI$(LIB_EXT) $(OBJECT)
	$(RANLIB) libPARI$(LIB_EXT)

' . qq{
pariinl.h: @inlines
	\$(PERL) -pe1 @inlines > pariinl.h

$converted

mp\$(OBJ_EXT): $mp pariinl.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) -o \$@ $mp

kernel\$(OBJ_EXT): $kernel1c pariinl.h
	\$(ASSCMD) \$(DEFINE) -o \$@ $kernel1c

kernel2\$(OBJ_EXT): $kernel2c pariinl.h
	\$(ASSCMD) \$(DEFINE) -o \$@ $kernel2c

mpinl\$(OBJ_EXT): $mpinl $mpinlh pariinl.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) -o \$@ $mpinl

} .
    join "\n", map { <<EOF } keys %cfiles;
$_\$(OBJ_EXT): $cfiles{$_} pariinl.h paricfg.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) -o \$@ $cfiles{$_}
EOF
}



sub MY::const_config
{
 my $self = shift;
 my $flags = $self->{'CCCDLFLAGS'}; # Tmp var needed with Perl4 !
 $flags =~ s/(-[fK]?\s*)pic\b/${1}PIC/; 
 $flags =~ s/-KPIC/-K PIC/;	# Apparently needed on Solaris...
 $self->{'CCCDLFLAGS'} = $flags;
 if ($^O eq 'MSWin32' && $Config{'ccflags'} =~ /-DPERL_OBJECT/)
  {
   $self->{'LDFLAGS'} =~ s/-(debug|pdb:\w+)\s+//g;
   $self->{'LDDLFLAGS'} =~ s/-(debug|pdb:\w+)\s+//g;
  }
 return $self->MM::const_config;
}

