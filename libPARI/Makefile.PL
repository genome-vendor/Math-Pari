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

%Versions = (
	     sun3	   => ["version68k.c",'sparc$(OBJ_EXT)'],
	     sparcv7	   => ["versionsparcv7.c",'sparcv7$(OBJ_EXT)'],
	     sparcv8micro  => ["versionsparcv8.c",'sparcv8micro$(SOL_UNDER)$(OBJ_EXT)'],
	     sparcv8super  => ["versionsparcv8.c",'sparcv8super$(SOL_UNDER)$(OBJ_EXT)'],
	     alpha	   => ["versionalpha.c",'alpha$(OBJ_EXT)'],
	     hppa	   => ["versionhp.c",'hppa$(OBJ_EXT)'],
	     i386	   => ["version386.c",'i386$(OBJ_EXT)'],
	     port	   => ["versionport.c"],
	    );

@Chips = keys %Versions;

@know_chip{@Chips} = (1) x @Chips;
$chipre = join '|', @Chips;
$chipre =~ s/\|port\b//;
$chipre =~ s/\bi386\b/[xi]\d86/;

if (not defined $machine and defined $common::machine) {
  $machine = $common::machine;
} elsif (not defined $machine and @ARGV and $ARGV[0] =~ /^machine=(.*)$/i) {
  $machine = $1; shift;
}

$machine = process_machine unless defined $machine;

if (not defined $machine 
    and $Config{myuname} 
    =~ /\b(sun3|sparcv7|sparcv8micro|sparcv8super|alpha|hppa|[ix]\d86)\b/) {
  $machine = $1;  
} elsif (not defined $machine) {
  chomp($machine = `uname -m`);
}
$machine =~ s/[ix]\d86/i386/ if defined $machine;

print "I think that your processor is of class `$machine'\n" 
  if defined $machine;

$sparc_cvt = $Config{osname} =~ /^(linux|nextstep)$/;

%inlines = (
	     sun3	   => ['none/level0.h','none/level1.h'],
	     sparcv7	   => ['none/asm0.h','none/level1.h'],
				# Following 2 guys should be fixed on linux/NS
	     sparcv8micro  => ( $sparc_cvt
				? ['none/asm0.h','none/level1.h']
				: ['sparcv8/level0.h','none/level1.h']),
	     sparcv8super  => ( $sparc_cvt
				? ['none/asm0.h','none/level1.h']
				: ['sparcv8/level0.h','none/level1.h']),
	     alpha	   => ['none/asm0.h','none/level1.h'],
	     hppa	   => ['none/asm0.h','none/level1.h'],
	     i386	   => ['ix86/level0.h','none/level1.h'],
	     port	   => ['none/level0.h','none/level1.h'],
	    );
$inlines = $inlines{$machine} or die "Unknown inlines for '$machine'";
@inlines = map "$main_paridir/src/kernel/$_", @$inlines;

%level0 = (
	     sun3	   => '',
	     sparcv7	   => '',
				# Following 2 guys should be fixed on linux/NS
	     sparcv8micro  => ($sparc_cvt 
			       ? ['sparc8/level0.s', 1,
				  'sparc8/level0_sparcv8_micro.s', 1] 
			       : ['sparc8/level0.s', 0, 
				  'sparc8/level0_sparcv8_micro.s', 0]),
	     sparcv8super  => '',
	     alpha	   => '',
	     hppa	   => ['none/level0.c', 0],
	     i386	   => ['ix86/l0asm.c', 1],
	     port	   => '',
	    );

%asmarch = (
	     sun3	   => 'm86k',
				# Following 2 guys should be fixed on linux/NS
	     sparcv8micro  => 'sparcv8',
	     sparcv8super  => 'sparcv8',
	     i386	   => 'ix86',
	     port	   => 'none',
	    );

$asmarch = $asmarch{$machine} || $machine;

$level0machine = $level0{$machine} || ["$asmarch/level0.s", 0];

$level0		= "$main_paridir/src/kernel/$level0machine->[0]";
$level0convert	= $level0machine->[1];
$level0c	= $level0convert ? 'level0.s' : $level0;
($level0h	= $level0) =~ s/\.c$/.h/;
$level0h	= '' unless -r $level0h;
($level0path	= $level0) =~ s/[^\/]*\.c$//;
$kernel1	= $kernel1c = $level0machine->[2];
$kernel1c	= 'kernel1.s'
  if $level0machine->[3];  

%mp = (
	     sun3	   => 'none/mp.c',
	     sparcv7	   => 'none/mp.c',
	     sparcv8micro  => 'none/mp.c',
	     sparcv8super  => 'none/mp.c',
	     alpha	   => 'none/mp.c',
	     hppa	   => 'none/mp.c',
	     i386	   => 'none/mp.c',
	     port	   => 'none/mp.c',
	    );
$mp = "$main_paridir/src/kernel/$mp{$machine}";

%mpinl = ( );
$mpinl_machine = $mpinl{$machine} || 'none/level1.c';
$mpinl = "$main_paridir/src/kernel/$mpinl_machine";
($mpinlh = $mpinl) =~ s/\.c$/.h/;

$libs="-lm";
$sol_under = ''; 

#@cfiles = ();
#opendir DIR, '.' or die "Cannot opendir . : $!";
#while (defined ($_ = readdir DIR)) {
#  push @cfiles, $_ if /\.c$/i 
#    and not /^(plot|i386|version|mpin|dummy)/;
#}
#closedir DIR or die "Cannot closedir . : $!";

@cfiles = <$main_paridir/src/*/*.c>;
@cfiles = grep !/\/(plot|i386|version|mpin|dummy)/, @cfiles;
%cfiles = map {m,/([^/.]*)\.c,i; ($1,$_)} @cfiles;


if ($Config{osname} eq 'solaris') {
  @sc_dirs = </opt/SUNWspro/SC*>;
  die "Cannot find SUNWspro dirs, needed for -lsunmath" unless @sc_dirs;
  $libs .= " -L$sc_dirs[-1]/lib -lsunmath";
  $sol_under = '_nounder';	# Need to filter
} elsif ($Config{osname} eq 'sunos') {
  # $machine = 'sparcv8super' if $Config{myuname} =~ /\bsun4m\b/o;
  $mycflags .= ' -DULONG_NOT_DEFINED';
}


unless (exists $know_chip{$machine}) {
  warn 
    "Do not know how make for chip `$machine', forcing portable interface\n";
  $machine = 'port';
}

# We remove optimize options, since Perl probably knows them
# and one can specify "OPTIMIZE=-O2 -m486" on the command line

if ($Config{cc} eq 'gcc') {
  #$mycflags="-O2 -g";		# We do not put -ansi to avoid versioncflags
  $mycflags="";		# We do not put -ansi to avoid versioncflags
} elsif ($machine eq 'hppa') {
  #$mycflags="-O -Aa -DHPPA";
  $mycflags="-Aa -DHPPA";
} elsif ($machine eq 'i386') {
  #$mycflags="-O2 -m486";
  $mycflags="";
} else {
  #$mycflags="-O";
  $mycflags="";
} 

$mycflags .= ' -Dshifts=pari_shifts'; # Conflict with libnsl
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

#@plotfiles = ('plotgnuplot.c', 'plotport.c');
$cfiles{plotgnuplot} = "$main_paridir/src/graph/plotgnuplot.c";
$cfiles{plotport} = "$main_paridir/src/graph/plotport.c";

#push @cfiles, $Versions{$machine}[0], @plotfiles, @portfiles;
#@obj_files = map {($a = $_) =~ s/\.c$/\$(OBJ_EXT)/; $a} @cfiles;
#push @obj_files, $Versions{$machine}[1] if defined $Versions{$machine}[1];

@obj_files = map { "$_\$(OBJ_EXT)" } keys %cfiles;
push @obj_files, 'kernel$(OBJ_EXT)', 'mp$(OBJ_EXT)', 'mpinl$(OBJ_EXT)', 
  ($kernel1c 
   ? 'kernel1$(OBJ_EXT)'
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

WriteMakefile(
    NAME	=> 'Math::PARI::libPARI',
    LINKTYPE	=> 'static',
    LIBS	=> $libs,
    OBJECT	=> join(' ', @obj_files), 
    macro	=> { 
		    MP_CFLAGS	 => $mp_cflags,
		    SOL_UNDER	 => $sol_under, 
		    ADD_AR_OPT	 => ($add_ar_flags || ''),
		    CPP => $Config{cppstdin},
		    CPPMINUS => $Config{cppminus},
		   },
    DEFINE	=> $mycflags,
    INC		=> "-I$main_paridir/src/headers -I$main_paridir/src/graph",
    C		=> \@cfiles,
    SKIP	=> [qw( distclean test dist makeaperl xs_o static)],
    clean	=> {'FILES' => 
		    'libPARI$(LIB_EXT) pariinl.h level0.s kernel1.s'},
);

sub MY::top_targets {
  '
all :: libPARI$(LIB_EXT)

static ::       libPARI$(LIB_EXT)

libPARI$(LIB_EXT): $(OBJECT)
	-$(RM_F) libPARI$(LIB_EXT)
	$(AR) $(ADD_AR_OPT) cr libPARI$(LIB_EXT) $(OBJECT)
	$(RANLIB) libPARI$(LIB_EXT)

' . qq{
pariinl.h: @inlines
	\$(PERL) -pe1 @inlines > pariinl.h

level0.s: $level0 $level0h
	\$(CPP) -I. \$(INC) -I$level0path \$(DEFINE) \$(CPPMINUS) < $level0 | perl -ne "s/%\s+/%/g; print unless /^\s*#/" > \$@

kernel1.s: $kernel1
	\$(CPP) -I. \$(INC) -I$level0path \$(DEFINE) \$(CPPMINUS) < $kernel1 | perl -ne "s/%\s+/%/g; print unless /^\s*#/" > \$@

mp\$(OBJ_EXT): $mp pariinl.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) -o \$@ $mp

kernel\$(OBJ_EXT): $level0c pariinl.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) -o \$@ $level0c

kernel1\$(OBJ_EXT): $kernel1c pariinl.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) -o \$@ $kernel1c

mpinl\$(OBJ_EXT): $mpinl $mpinlh pariinl.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) -o \$@ $mpinl

} . 
    join "\n", map { <<EOF } keys %cfiles;
$_\$(OBJ_EXT): $cfiles{$_} pariinl.h paricfg.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) -o \$@ $cfiles{$_}
EOF
}


