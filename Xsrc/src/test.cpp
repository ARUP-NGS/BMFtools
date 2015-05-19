#include "kseq.hpp"
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int main(int argc, char * argv[]) 
{
    if(argc != 3) {
        std::cout<< argv[0]<<" "<<" input_file compression"<<std::endl;
        std::cout<<"e.g. "<< argv[0]<<" file.fa gzip"<<std::endl;
        std::cout<<"e.g. "<< argv[0]<<" file.fa bzip2"<<std::endl;
        std::cout<<"e.g. "<< argv[0]<<" file.fa none"<<std::endl;
    }
    kseq seq;
    int l = 0;
    int c = 0;
    switch (argv[2][0]) 
    {
        case 'g':
          {
          gzFile fp = gzopen(argv[1], "r");
          FunctorZlib gzr;
          kstream<gzFile, FunctorZlib> ks(fp, gzr);
          while((l = ks.read(seq)) >= 0) {
            std::cout << seq.name << std::endl;
            std::cout << seq.seq << std::endl;
            c++;
          }
          gzclose(fp);
          break;
          }
        case 'b':
          {
            BZFILE * fp = BZ2_bzopen(argv[1], "r");
            FunctorBZlib2 bzr;
            kstream<BZFILE*, FunctorBZlib2>ks(fp, bzr);
            while((l = ks.read(seq)) >= 0) {
              std::cout << seq.name << std::endl;
              std::cout << seq.seq << std::endl;
              c++;
            }
            BZ2_bzclose(fp);
            break;
          }
        default:
          {
            int fp = open(argv[1], O_RDONLY);
            FunctorRead r;
            kstream<int, FunctorRead> ks(fp, r);

            while((l = ks.read(seq)) >= 0) {
              std::cout << seq.name << std::endl;
              std::cout << seq.seq << std::endl;
              c++;
            }
            close(fp);
          }
    }



    std::cout << c << std::endl;

    return 0;
}
