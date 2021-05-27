#!/usr/bin/env perl
# gffggb.perl
# cat dmel_3L_172k.thmap2.gff | perl -pe's/ID=\w+;//' | sort -k3,3 -k9,9 | env s=src1 perl gfgb.perl \
# | sort -k3,3n > dmel_3L_172k.thmap2.ggb

BEGIN{ $st=$ENV{'s'} || "rna"; }
$dropsub=$ENV{nosub}||0;

while(<>){
next unless(/^\w/); @v=split"\t"; 
($g)=m/Parent=([^;\s]+)/; ($r,$t,$b,$e,$o)=@v[0,2,3,4,6];
$g="g".++$gn unless($g); 
$g =~ s/(\d)t\d+$/$1/  if $t eq "intron";  # drop alttr, same intron
next if($t eq "intron" and $didi{$g.$b.$e.$o}++);
$g.= ".". ++$in if $t eq "intron"; 
if($lg and ($lg ne $g or $lt ne $t)){ put(); } 
push(@be,"$b-$e"); $lo=$o unless($o eq "."); $lt=$t; $lg=$g; $lr=$r; 
}

END{ put(); ggbhead(); puto();} 

sub revo() { map{ s/(\d+).(\d+)/$2-$1/ }@be; @be=sort{$b<=>$a}@be; } 

sub put{ revo() if ($lo eq "-"); $gtag="$st$lt" unless($gtag); $vo= "$st$lt\t$lg\t". join(",",@be)."\n"; 
push @vo, $vo;  @be=(); $lo="."; }

sub puto{ print "\nreference = $lr\n\n" unless($didr{$lr}++); 
unless($dropsub) { foreach $vo (@vo) { print $vo; } }
else {
foreach (@vo) { ($t,$d,$x)=split; $xd{$d}=$x; @x=split /[,-]/,$x;
($b,$e)= @x[0,-1]; ($b,$e)=($e,$b) if($b>$e); $w=1+$e-$b; $xw{$d}=$w; }
@d=sort{$xw{$b}<=>$xw{$a}} keys %xw; foreach $d (@d) { $x=$xd{$d};
($s=$x)=~s/^(\d+)//; $xb=$1; $s=~s/(\d+)$//; $xe=$1; if(($i=index($xs,$s))>=0) { }
else { $xs.="$d:$x;" } } 
@xs=split";",$xs; foreach $xs (@xs) { ($d,$x)=split":",$xs; 
print join("\t",$t,$d,$x)."\n" ; } 
} }


## add filter for subset genes:
# perl -ne'($t,$d,$x)=split; unless(/^v9c/) { print; next; }  $xd{$d}=$x; @x=split /[,-]/,$x;
# ($b,$e)= @x[0,-1]; ($b,$e)=($e,$b) if($b>$e); $w=1+$e-$b; $xw{$d}=$w; 
# END{ @d=sort{$xw{$b}<=>$xw{$a}} keys %xw; foreach $d (@d) { $x=$xd{$d}; 
# ($s=$x)=~s/^(\d+)//; $xb=$1; $s=~s/(\d+)$//; $xe=$1; if(($i=index($xs,$s))>=0) { } 
# else { $xs.="$d:$x;" } } @xs=split";",$xs; foreach $xs (@xs) {
#  ($d,$x)=split":",$xs; print join("\t",$t,$d,$x)."\n" ; } }' 


sub colr{ my @co=qw(coral darkslateblue green lightpink mediumslateblue paleturquoise sienna
black cornflowerblue darkslategray greenyellow lightsalmon mediumspringgreen palevioletred
silver aliceblue cornsilk darkturquoise honeydew lightseagreen mediumturquoise papayawhip
skyblue crimson darkviolet hotpink lightskyblue mediumvioletred peachpuff slateblue aqua
cyan deeppink indianred lightslategray midnightblue peru slategray aquamarine darkblue
deepskyblue indigo lightsteelblue mintcream pink snow azure darkcyan dimgray ivory
lightyellow mistyrose plum springgreen beige darkgoldenrod dodgerblue khaki lime moccasin
powderblue steelblue bisque darkgray firebrick lavender limegreen purple tan
blanchedalmond darkgreen lavenderblush linen navy red teal blue darkkhaki forestgreen
lawngreen magenta oldlace rosybrown thistle blueviolet darkmagenta fuchsia lemonchiffon
maroon olive royalblue tomato brown darkolivegreen gainsboro lightblue mediumaquamarine
olivedrab saddlebrown turquoise burlywood darkorange lightcoral mediumblue orange salmon
violet cadetblue darkorchid gold lightcyan mediumorchid orangered sandybrown wheat
chartreuse darkred goldenrod lightgoldenrodyellow mediumpurple orchid seagreen chocolate
darksalmon gray lightgreen mediumseagreen palegoldenrod seashell yellow coral darkseagreen
green lightgrey mediumslateblue palegreen sienna yellowgreen); return $co[ int(rand(@co))]; }

## ggb header
sub ggbhead { 
$stx=$gtag; unless($stx) { $stx=$ENV{'s'}||"rna"; $stx.="exon"; }
$colr=colr();
print <<"EOT";
[$stx]
glyph = transcript2
strand = 1
height = 6
bgcolor = $colr
bump = 1

EOT
}

