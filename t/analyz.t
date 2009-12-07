$file = __FILE__;
$file =~ m|([^/.]*)\.t$| or die;
@ARGV = "../src/test/32/$1";
do 'test_eng/Testout.t';
die if $@;

