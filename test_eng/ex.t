$file = __FILE__;
$file =~ m|([^/.]*)\.t$| or die;
@ARGV = "CHANGE_ME/src/test/32/$1";
do 'test_eng/Testout.pm';
die if $@;
