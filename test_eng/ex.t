$file = __FILE__;
$file =~ m|^(.*)[\\/]([^\\/.]*)\.t$|s or die;
$dir = $1;
$name = $2;
$dir1 = "CHANGE_ME";
$dir1 = "$dir/../$dir1" unless $dir1 =~ m|^([a-z]:)?[\\/]|i;
@ARGV = "$dir1/src/test/32/$name";
do 'test_eng/Testout.pm';
die if $@;
