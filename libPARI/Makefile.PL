use ExtUtils::MakeMaker;
use Config '%Config';

#BEGIN {unshift @INC, 'utils'}
use Math::PariBuild;

ExtUtils::MakeMaker::WriteEmptyMakefile(), return 1 if $common::parilib;

*MY::top_targets  = \&xMY::top_targets;		# will break WriteEmptyMakefile
*MY::const_config = \&xMY::const_config;

$main_paridir = $common::main_paridir;
die "Please start Makefile.PL at toplevel, not in libPARI"
  unless $main_paridir;

if ($main_paridir =~ m,^([a-z]:)?[/\\],i) {
  $main_paridir = $common::main_paridir;
} else {
  $main_paridir = "../$common::main_paridir";
}

$pari_version = $common::pari_version;
$machine = $common::machine unless defined $machine;

if (defined $machine) {
  print "...Via command-line: processor family `$machine'\n";
} else {
  $machine = find_machine_architecture;
}
if ( $machine eq 'port'			# Set by find_machine_architecture()
     and defined($Config{longsize}) and $Config{longsize} == 8 ) {
  #$mycflags .= " -DLONG_IS_64BIT";	# Moved to BuildPari...
}

# Choose which assembler files to use
$asmarch = choose_and_report_assembler($machine);
# $machine used as a coarse-grain factor below, $asmarch as a fine-grain
$machine = 'port' if $asmarch eq 'none';

# This part is based on analysing Makefile.SH in */kernel/* subdirectories.

if ($machine =~ /sparc|alpha|hppa/) {
  $asscmd = '$(AS) $(ASFLAGS)';	#  $(CCCDLFLAGS)?
} else {
  $asscmd = '$(CCCMD) $(CCCDLFLAGS)';
}

$Using_gnu_as = is_gnu_as();
print STDERR "...Assembler is GNU assembler\n" if $Using_gnu_as;

$asmcpp_define = '';
if ($Using_gnu_as and $^O eq 'solaris' and $asmarch =~ /^sparc/) {
  $asmcpp_define = '-DNO_UNDERSCORE';
  $asmcpp_define .= ' -D__GNUC__' if $Config{gcc};
}

# Which files to compile and how to build levels1/2 kernel
$kernels = kernel_files($asmarch, $pari_version, $Using_gnu_as);
kernel_fill_data($kernels, \%kernel);

# remaining kernel C files
$kernel_dir  = '$(PARI_DIR)/src/kernel';
$mp	     = "$kernel_dir/none/mp.c";
$mpinl	     = "$kernel_dir/none/level1.c"; # Non-inlined versions of inlines
$mpinlh	     = "$kernel_dir/none/level1.h"; # Put as a dependence


@cfiles =			# Remaining C files
  grep !m(/ ( plot (?!port|gnuplot)
	    | ix86 | version | mpin | dummy
	    | gp_(rl|init) | kerntest | whatnow | gp\. )
	 )x,	<$main_paridir/src/*/*.c>, # kernel C files deeper in the tree
		<$main_paridir/src/systems/$^O/*.c>;
@cat_cfiles =			# Remaining kernel C files
  grep !m(/ ( level \d+ | mp ) )x,  <$main_paridir/src/kernel/none/*.c>;


@cfiles = grep ! m(/test/[^/]+\.c$), @cfiles;
map s/^\Q$main_paridir/\$(PARI_DIR)/, @cfiles;
%cfiles = map {m,/([^/.]*)\.c$,i; ($1,$_)} @cfiles; # By basename

$noexp2 = '';			# Not used any more - but would be nice!
$libs="-lm";
if ($Config{osname} eq 'solaris') {
#  @sc_dirs = '/opt/SUNWspro/lib' if -d '/opt/SUNWspro/lib';
#  @sc_dirs = </opt/SUNWspro/SC*/lib> unless @sc_dirs;
  #warn "Cannot find SUNWspro dirs, needed for -lsunmath, using NOEXP2.\n"
  #  unless @sc_dirs;
#  if (@sc_dirs) {
#    $libs .= " -L$sc_dirs[-1] -lsunmath";
#  } else {
    $noexp2 = 1;
#  }
} elsif ($Config{osname} eq 'sunos') {
  # $machine = 'sparcv8_super' if $Config{myuname} =~ /\bsun4m\b/o;
  $mycflags .= ' -DULONG_NOT_DEFINED';
}

# We remove optimize options, since Perl probably knows them
# and one can specify "OPTIMIZE=-O2 -m486" on the command line

my $extra_inc = extra_includes($paridir);

my $inc = '';
if ($Config{cc} eq 'gcc') {
  # We do not put -ansi to avoid versioncflags
  $mycflags="";			#$mycflags="-O2 -g";
  if ($Config{osname} eq 'os2' and not $extra_inc) {
    # 2.1.3 does not include dlctl.h
    $mycflags .= ' -DRTLD_LAZY=0 -DRTLD_GLOBAL=0';
    $inc .= ' -I$(PERL_INC)';
  }
} elsif ($machine eq 'hppa') {
  $mycflags="-Aa -DHPPA";	#$mycflags="-O -Aa -DHPPA";
} elsif ($machine eq 'ix86') {
  $mycflags="";			#$mycflags="-O2 -m486";
} else {
  $mycflags="";			#$mycflags="-O";
}

# These are not needed with newer versions of PARI
# $mycflags .= ' -Dshifts=pari_shifts'; # Conflict with libnsl

$mycflags .= ' -Derr=pari_err'; # On linux it can get a wrong dynamic loading

if ($machine eq 'alpha') {
  #$mycflags .= " -DLONG_IS_64BIT";	# Moved to BuildPari...
} elsif ($Config{osname} eq 'solaris') {
  $mycflags .= " -DSOLARIS";
} elsif ($Config{osname} eq 'os2') {
  $noexp2 = 1;
  $mycflags .= ' -DMALLOC_PROCS -DHAS_STRICMP';
  $add_ar_flags = '-p64';	# Need big library. 32 is OK without debugging
} elsif ($Config{osname} eq 'linux') {
  $noexp2 = 1;
}

if ($Config{gccversion}
    and not $common::do_configure
    and not ($asmarch eq 'alpha'
	     and $Config{gccversion} !~ /^(3\.)|(2\.95\.3.*)/)) {
  $mycflags .= " -DASMINLINE";
}
$mycflags .= ' -DGCC_INLINE' if $Config{gccversion};

@obj_files = map { "$_\$(OBJ_EXT)" } keys(%cfiles), qw(kernel mp mpinl);
push @obj_files, 'kernel2$(OBJ_EXT)' if $kernel{converted2};

$mycflags .= " -DDYNAMIC_PLOTTING";
# OMF build needs no underscores:
$mycflags .= " -D__NO_AOUT" if $^O eq 'os2' and !$OS2::is_aout;

# Fix for broken gcc 2.95 with -ffast-math
my @opt;
my $opt = $Config{optimize};

if (defined $Config{gccversion} and $Config{gccversion} =~ /^2\.9[46]/
    and $opt =~ s/-ffast-math// and not "@ARGV" =~ /\bOPTIMIZE=/) {
  @opt = ( OPTIMIZE => $opt );
}

my $Using_ms_vc = ($^O =~ /win32/i and $Config{cc} =~ /cl/i);
my $Using_Borland = ($^O =~ /win32/i and $Config{cc} =~ /\bbcc/i);

WriteMakefile(
    NAME	=> 'Math::PARI::libPARI',
    LINKTYPE	=> 'static',
    LIBS	=> $libs,
    OBJECT	=> join(' ', @obj_files),
    @opt,
    macro	=> {
		    PARI_DIR	  => $main_paridir,
		    ADD_AR_OPT	  => ($add_ar_flags || ''),
		    CPP		  => $Config{cpp},
		    # CPPMINUS	  => $Config{cppminus},
		    ASSCMD	  => $asscmd,
				# Extra defines when CPPing assembler files:
		    ASSCPPDEFINE  => $asmcpp_define,
				# Use this when assembling:
		    ASFLAGS	  => assembler_flags($machine, $Using_gnu_as),
				# defines used when assembling:
		    ASSDEFINE	  => ($Using_gnu_as ? '' : '$(DEFINE)'),
		    MY_CC_PRE_TARGET => ($Using_ms_vc
					 ? '-Fo'
					 : ($Using_Borland
					    ? '-o': '-o $(MY_EMPTY_STR)')),
		    MY_AR_PRE_TARGET => ($Using_ms_vc
					 ? '-out:'
					 : ($Using_Borland
					    ? '' : 'cr $(MY_EMPTY_STR)')),
		    MY_AR_OBJECT => ($Using_Borland
				     ? '$(OBJECT:^"+")'
				      : '$(OBJECT)'),
		    MY_EMPTY_STR  => '',
		   },
    DEFINE	=> $mycflags,
    INC		=> $extra_inc . ' -I$(PARI_DIR)/src/headers -I$(PARI_DIR)/src/graph -I.'
			. $inc,
    C		=> \@cfiles,
    SKIP	=> [qw( distclean test dist makeaperl xs_o static)],
    clean	=> {FILES =>
		    'libPARI$(LIB_EXT) pariinl.h kernel1.s kernel2.s paricfg.h mp.c kernel1.c'},
);

# Create preprocessed copies of assembler files (GNU as does not run cpp);
# Create .h file with inlining rules and an object file mpinl with
# non-inlined versions of normally inlined functions.
# Assemble assembler files.
# Compile C files from a *different* directory.
# Create a library;
sub xMY::top_targets {
  my $converted = '';
  my @inlines = inline_headers($asmarch, $pari_version);
  $converted = <<EOC if $kernel{convert};
$kernel{converted1}: $kernel{file1} $kernel{header1}
	\$(CPP) -I. \$(INC) -I$kernel{dir1} \$(DEFINE) \$(ASSCPPDEFINE) $kernel{file1} | perl -ne "s/%\\s+/%/g; print unless /^\\s*#/" > \$@

kernel2.s: $kernel{file2}
	\$(CPP) -I. \$(INC) -I$kernel{dir1} \$(DEFINE) \$(ASSCPPDEFINE) $kernel{file2} | perl -ne "s/%\\s+/%/g; print unless /^\\s*#/" > \$@

EOC
  my $create_mp = '';
  if ($pari_version >= 2002005) {
    my $mp_pre = $mp;
    $mp = 'mp.c';
    my $kern = "$kernel_dir/none";
    my @dep			# Copied from src/kernel/none/MakeLVL1.SH 2.2.10
      = "$kern/mp.c $kern/cmp.c $kern/gcdll.c $kern/ratlift.c $kern/gcd.c $kern/invmod.c $kern/mp_indep.c $kern/add.c";
    $create_mp = <<EOS;
mp.c: @dep
	\$(PERL) -pe1 @dep > \$@
EOS
  }
  return '
all :: libPARI$(LIB_EXT)

static ::       libPARI$(LIB_EXT)

test:

libPARI$(LIB_EXT): $(OBJECT)
	-$(RM_F) libPARI$(LIB_EXT)
	$(AR) $(ADD_AR_OPT) $(MY_AR_PRE_TARGET)libPARI$(LIB_EXT) $(MY_AR_OBJECT)
	$(RANLIB) libPARI$(LIB_EXT)

' . qq{
$create_mp
pariinl.h: @inlines
	\$(PERL) -pe1 @inlines > pariinl.h

$converted

mp\$(OBJ_EXT): $mp pariinl.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) \$(MY_CC_PRE_TARGET)\$@ $mp

kernel\$(OBJ_EXT): $kernel{converted1} pariinl.h
	\$(ASSCMD) \$(ASSDEFINE) \$(MY_CC_PRE_TARGET)\$@ $kernel{converted1}

kernel2\$(OBJ_EXT): $kernel{converted2} pariinl.h
	\$(ASSCMD) \$(ASSDEFINE) \$(MY_CC_PRE_TARGET)\$@ $kernel{converted2}

mpinl\$(OBJ_EXT): $mpinl $mpinlh pariinl.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) \$(MY_CC_PRE_TARGET)\$@ $mpinl

} .
    join "\n", map { <<EOF } sort keys %cfiles;
$_\$(OBJ_EXT): $cfiles{$_} pariinl.h paricfg.h
	\$(CCCMD) \$(CCCDLFLAGS) \$(DEFINE) \$(MY_CC_PRE_TARGET)\$@ $cfiles{$_}
EOF
}


# Rewrite pic option to PIC in CCCDLFLAGS,
# Why manipulate LD* when no ld is done here (static build only?)
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

