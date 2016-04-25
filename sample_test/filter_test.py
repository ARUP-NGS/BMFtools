#!/usr/bin/env python

if __name__ == "__main__":

    import unittest
    from subprocess import check_output

    bampath = "FINAL_S376x3_Unpurified_GCATAACG-GCATAACG_L001_R1_001.R1.rsq.bam"
    class Test(unittest.TestCase):
        def test_empty_filter(self):
            pre_count = check_output("samtools view -cF2304 %s" % bampath, shell=True)
            post_count = check_output("bmftools_db filter -l0 %s - | samtools view -cF2304" % bampath, shell=True)
            self.assertEqual(pre_count, post_count)

        def test_mq_filter(self):
            pre_count = check_output("samtools view -cF2304 -q30 %s" % bampath, shell=True)
            post_count = check_output("bmftools_db filter -m30 -l0 %s - | samtools view -cF2304" % bampath, shell=True)
            self.assertEqual(pre_count, post_count)

        def test_fm_filter(self):
            pre_count = check_output("samtools view -F2304 %s | grep -Pvc 'FM:i:1\t'" % bampath, shell=True)
            post_count = check_output("bmftools_db filter -F2304 -s2 -l0 %s - | samtools view -c" % bampath, shell=True)
            self.assertEqual(pre_count, post_count)

        def test_flag_filter(self):
            for flag in [4, 17, 2304, 1173]:
                pre_count = check_output("samtools view -cf%i %s" % (flag, bampath), shell=True)
                post_count = check_output("bmftools_db filter -f%i -l0 %s - | samtools view -c" % (flag, bampath), shell=True)
                self.assertEqual(pre_count, post_count)
                pre_count = check_output("samtools view -cF%i %s" % (flag, bampath), shell=True)
                post_count = check_output("bmftools_db filter -F%i -l0 %s - | samtools view -c" % (flag, bampath), shell=True)
                self.assertEqual(pre_count, post_count)

        def test_split_filter(self):
            pre_count = int(check_output("samtools view -c %s" % bampath, shell=True).strip())
            post_count = int(check_output("bmftools_db filter -F2304 -rtmpwtf.bam -s3 %s - | samtools view -c" % bampath, shell=True).strip())
            post_count += int(check_output("samtools view -c tmpwtf.bam", shell=True).strip())
            check_output("rm tmpwtf.bam", shell=True)
            self.assertEqual(pre_count, post_count)

    unittest.main()
