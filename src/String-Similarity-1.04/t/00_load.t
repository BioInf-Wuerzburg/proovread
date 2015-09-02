BEGIN { $| = 1; print "1..10\n"; }
END {print "not ok 1\n" unless $loaded;}
use String::Similarity;
$loaded = 1;
print "ok 1\n";

print similarity("this should be the same", "this should be the same") == 1 ? "" : "not ", "ok 2\n";
my $s = similarity("this should be same the", "this should be the same");
print $s > 0.825 && $s < 0.827 ? "" : "not ", "ok 3\n";
print similarity("A", "B") == 0 ? "" : "not ", "ok 4\n";
print similarity("\x{456}", "\x{455}") == 0 ? "" : "not ", "ok 5\n";

require utf8;

$a = chr(169);
utf8::upgrade($b = $a);
print similarity($a, $b) == 1 ? "" : "not ", "ok 6\n";

$b = "\x{0040}";
utf8::downgrade($a = $b);
print similarity($a, $b) == 1 ? "" : "not ", "ok 7\n";

utf8::upgrade($a = chr(169));
use Encode qw( encode );
$b = encode "utf-8", $a;
print similarity($a, $b) < 1 ? "" : "not ", "ok 8\n";

$a = []; $b = [];
$s1 = similarity($a, $b);
$s2 = similarity($a, "$b");
print "not " unless abs($s1-$s2) < 0.001;
print "ok 9\n";

$s = similarity("", undef);
print "not " unless $s == 1;
print "ok 10\n";
