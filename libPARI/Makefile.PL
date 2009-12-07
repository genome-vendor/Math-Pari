use ExtUtils::MakeMaker;
use Config '%Config';
# Based on Makemakefile shell script by Bruno Haible 10.12.1992,
# modified by Babecool 23.7.1993
# Created by Ilya Zakharevich <ilya@math.ohio-state.edu>

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

$solaris='';

$libs="-lm";
$sol_under = ''; 
@plotfiles = ('plotgnuplot.c', 'plotport.c');

@cfiles = ();
opendir DIR, '.' or die "Cannot opendir . : $!";
while (defined ($_ = readdir DIR)) {
  push @cfiles, $_ if /\.c$/i 
    and not /^(plot|i386|version|mpin|dummy)/;
}
closedir DIR or die "Cannot closedir . : $!";

if ($Config{osname} eq 'solaris') {
  @sc_dirs = </opt/SUNWspro/SC*>;
  die "Cannot find SUNWspro dirs, needed for -lsunmath" unless @sc_dirs;
  $libs .= " -L$sc_dirs[-1]/lib -lsunmath";
  $gccversioncflags=1;
  $solaris='-DSOLARIS';
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

if ($machine eq 'alpha') {
  $mycflags .= " -DLONG_IS_64BIT";
} elsif ($Config{osname} eq 'solaris') {
  $mycflags .= " -DSOLARIS";
} elsif ($Config{osname} eq 'os2') {
  $mycflags .= ' -DNOEXP2 -DMALLOC_PROCS -DHAS_STRICMP';
  $add_ar_flags = '-p64';	# Need big library. 32 is OK without debugging
} elsif ($Config{osname} eq 'linux') {
  $mycflags .= ' -DNOEXP2';
}

$mp_cflags = '';
if ($machine eq 'sun3') {
  $mp_cflags = '\$(CPlPl) \$(NOASM)'; # NOASM seems to be empty
} elsif ($machine eq 'port') {
  push @portfiles, 'mp.c';
  $mycflags .= " -D__HAS_NO_ASM__";
}

push @cfiles, $Versions{$machine}[0], @plotfiles, @portfiles;
@obj_files = map {($a = $_) =~ s/\.c$/\$(OBJ_EXT)/; $a} @cfiles;
push @obj_files, $Versions{$machine}[1] if defined $Versions{$machine}[1];

unless (-f 'genpari.h') {
  $genpari = ($machine eq 'sun3' ? 'genpari68k.h' : 'genpariother.h');
  eval {link $genpari, 'genpari.h'} or do {
    require File::Copy; 
    File::Copy::copy ($genpari, 'genpari.h') 
      or die "Cannot copy $genpari to genpari.h";
  }
}

$mycflags .= " -DDYNAMIC_PLOTTING";

WriteMakefile(
    'NAME'	=> 'Math::PARI::libPARI',
    'LINKTYPE'	=> 'static',
    'LIBS'	=> $libs,
    'OBJECT'	=> join(' ', @obj_files), 
    macro	=> { 
		    MP_CFLAGS	 => $mp_cflags,
		    SOL_UNDER	 => $sol_under, 
		    ADD_AR_OPT	 => ($add_ar_flags || ''),
		   },
    'DEFINE'	=> $mycflags,
    'C'		=> \@cfiles,
    'SKIP'	=> [qw( distclean test dist makeaperl xs_o static)],
    'clean'	=> {'FILES' => 'libPARI$(LIB_EXT)'}
);

sub MY::top_targets {
  '
all :: libPARI$(LIB_EXT)

static ::       libPARI$(LIB_EXT)

libPARI$(LIB_EXT): $(OBJECT)
	$(AR) $(ADD_AR_OPT) cr libPARI$(LIB_EXT) $(OBJECT)
	$(RANLIB) libPARI$(LIB_EXT)

sparcv7$(OBJ_EXT):	sparcv7.s
	as -P -o sparcv7$(OBJ_EXT) sparcv7.s

sparcv8micro$(SOL_UNDER)$(OBJ_EXT):	sparcv8micro$(SOL_UNDER).s
	as -P -o sparcv8micro$(SOL_UNDER)$(OBJ_EXT) sparcv8micro$(SOL_UNDER).s

sparcv8super$(SOL_UNDER)$(OBJ_EXT):	sparcv8super$(SOL_UNDER).s
	as -P -o sparcv8super$(SOL_UNDER)$(OBJ_EXT) sparcv8super$(SOL_UNDER).s

sparcv8micro_nounder.s: sparcv8micro.s
	sed s/_//g $? > $@

sparcv8super_nounder.s: sparcv8super.s
	sed s/_//g $? > $@

alpha$(OBJ_EXT):	alpha.s
	as -O1 -o alpha$(OBJ_EXT) alpha.s

hppa$(OBJ_EXT):	hppa.s
	cc $(CCFLAGS) -c -o hppa$(OBJ_EXT) hppa.s

i386$(OBJ_EXT):	i386.s
	$(CC) $(CCFLAGS) -c i386.s

i386.s:	i386.c
	$(CC) $(CCFLAGS) -E i386.c > i386.s

' . ( $machine eq 'sun3' ?
'mp$(OBJ_EXT):	mp.s
	$(CC) $(CCFLAGS) $(MP_CFLAGS) -c mp.s

' : '');
}


