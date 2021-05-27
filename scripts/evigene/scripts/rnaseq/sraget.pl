#!/usr/bin/env perl
# mirror fetch from ncbi sra_result.csv table listing ftp address of data
# # $mr ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX151/SRX151669

## fail unless `which lftp`; alt use wget ?
$mr='lftp -c mirror ';
while(<>) {
  chomp; @v= map{ s/^"//; s/"$//; $_; } split",";
  $ok=0; ($u)= grep /ftp:/, @v;
  if($u) { 
    if(0) { $ok= system("$mr $u ".'&'); }
    else {
    $pid=fork();
    if ($pid) { push @pid, $pid; }
    else { $ok= system('lftp', '-c', 'mirror', $u); } 
    }
   }
  warn "#fork=$ok,$pid $mr $u\n"; 
}

wait;
# waitpid @pid;

