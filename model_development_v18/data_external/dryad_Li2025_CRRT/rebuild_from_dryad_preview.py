import io, sys, csv, os
sys.stdout.reconfigure(encoding="utf-8")

CAZ = """subjectID\tTime_h\tcaz_pre\tcaz_post
1\t0\t82.5\t58.5
1\t2\t178\t123
1\t3\t148\t110
1\t4\t128\t104
1\t6\t95\t74.5
1\t8\t59.5\t48
2\t0\t63.5\t38
2\t2\t211\t138
2\t3\t154\t116
2\t4\t131\t83.5
2\t6\t93.5\t64.5
2\t8\t52\t36.5
3\t0\t75.5\t48.5
3\t2\t161\t102
3\t3\t134\t86
3\t4\t116\t83
3\t5\t90.5\t69
3\t8\t72\t54
4\t0\t59\t53.5
4\t2\t140\t135
4\t3\t107\t104
4\t4\t83.5\t84.5
4\t6\t80\t74
4\t8\t62\t61.5
5\t0\t41\t31
5\t1\t107\t76.5
5\t2\t92\t69
5\t3\t80\t58
5\t4\t75\t54.5
5\t6\t58\t40
5\t8\t44.5\t31
6\t0\t89\t58
6\t2\t198\t125
6\t3\t145\t102
6\t4\t132\t75.5
6\t6\t98.5\t63
6\t8\t72\t47
7\t0\t56.5\t52.5
7\t2\t92.5\t93.5
7\t3\t69.5\t64
7\t4\t57.5\t56
7\t6\t44.5\t43.5
7\t8\t37.5\t33
8\t0\t41.5\t30
8\t1\t116\t69.5
8\t2\t77.5\t48
8\t3\t80.5\t53.5
8\t4\t83\t52.5
8\t6\t66.5\t45.5
8\t8\t52.5\t33
9\t0\t58.5\t58
9\t3\t109\t98
9\t4\t87\t84.5
9\t6\t67.5\t61
9\t8\t54\t52.5
10\t0\t65.5\t48.5
10\t2\t164\t109
10\t3\t132\t90
10\t4\t109\t76.5
10\t6\t86\t60.5
10\t8\t64\t49.5
11\t0\t47\t30
11\t2\t153\t105
11\t3\t97\t75
11\t4\t84\t58.5
11\t6\t70\t53.5
11\t8\t60\t57.5
12\t0\t59\t61
12\t3\t100\t97
12\t4\t89.5\t87.5
12\t6\t69.5\t71.5
12\t8\t61.5\t60
13\t0\t49\t36.5
13\t3\t121\t92.5
13\t4\t97.5\t69.5
13\t6\t77.5\t58
13\t8\t61.5\t44
14\t0\t59\t45.5
14\t3\t112\t80.5
14\t4\t98.5\t67
14\t6\t80.5\t53
14\t8\t58\t38.5
15\t0\t83\t53
15\t3\t137\t101
15\t4\t112\t78
15\t6\t96\t57
15\t8\t65\t48
16\t0\t111\t98.5
16\t3\t186\t145
16\t4\t154\t123
16\t6\t120\t93
16\t8\t105\t85
17\t0\t67.5\t57
17\t3\t156\t132
17\t4\t134\t98
17\t6\t112\t87
17\t8\t85.5\t76.5
18\t0\t96.5\t73
18\t3\t175\t138
18\t4\t147\t116
18\t6\t114\t92
18\t8\t89.5\t72
19\t0\t60.5\t46
19\t1\t172\t124
19\t2\t158\t112
19\t3\t133\t93.5
19\t4\t119\t85.5
19\t6\t93\t67.5
19\t8\t73\t53
20\t0\t62\t61.5
20\t3\t139\t137
20\t4\t117\t118
20\t6\t87\t90
20\t8\t78\t78
21\t0\t74\t71
21\t3\t139\t138
21\t4\t131\t129
21\t6\t107\t105
21\t8\t76\t77"""

AVI = """subjectID\tTime_h\tavi_pre\tavi_post
1\t0\t16\t11
1\t2\t32.5\t22
1\t3\t27\t21.5
1\t4\t23.5\t21
1\t6\t18.5\t13.5
1\t8\t11.5\t9.5
2\t0\t14\t8.5
2\t2\t38\t26.5
2\t3\t29\t21.5
2\t4\t24.5\t16
2\t6\t19.5\t13
2\t8\t11\t8
3\t0\t13.5\t9
3\t2\t31\t18.5
3\t3\t25.5\t15.5
3\t4\t22.5\t16
3\t5\t18\t12.5
3\t8\t15\t10
4\t0\t14\t12.5
4\t2\t27.5\t27.5
4\t3\t22.5\t21.5
4\t4\t18\t17.5
4\t6\t17\t16.5
4\t8\t14.5\t13.5
5\t0\t8\t5.5
5\t1\t21.5\t13
5\t2\t18\t12.5
5\t3\t16\t12.5
5\t4\t15.5\t10.5
5\t6\t11.5\t7.5
5\t8\t9.5\t6.5
6\t0\t19.5\t11.5
6\t2\t37\t21.5
6\t3\t29\t19.5
6\t4\t26\t15
6\t6\t20\t13
6\t8\t15.5\t10
7\t0\t11.5\t9.5
7\t2\t28.5\t27
7\t3\t23.5\t21
7\t4\t20\t18.5
7\t6\t15\t15
7\t8\t14.5\t14.5
8\t0\t8\t5
8\t1\t21\t12.5
8\t2\t14.5\t9
8\t3\t15\t10
8\t4\t15\t9.5
8\t6\t13\t8
8\t8\t9.5\t5.5
9\t0\t13.5\t12.5
9\t3\t23.5\t20
9\t4\t18.5\t18.5
9\t6\t16.5\t14.5
9\t8\t12.5\t11.5
10\t0\t13\t9
10\t2\t29.5\t19.5
10\t3\t23.5\t16
10\t4\t18.5\t15
10\t6\t15\t11
10\t8\t12\t9
11\t0\t8.5\t5.5
11\t2\t29.5\t20.5
11\t3\t22.5\t13.5
11\t4\t20\t11
11\t6\t16.5\t9.5
11\t8\t13.5\t11.5
12\t0\t14\t14
12\t3\t23\t21.5
12\t4\t20\t20
12\t6\t16.5\t16
12\t8\t14\t13.5
13\t0\t13\t9
13\t3\t27\t21.5
13\t4\t23.5\t16.5
13\t6\t18.5\t14
13\t8\t14.5\t10.5
14\t0\t15.5\t11.5
14\t3\t29\t20.5
14\t4\t26.5\t18
14\t6\t21\t14
14\t8\t16\t10.5
15\t0\t20\t12.5
15\t3\t34\t23
15\t4\t26.5\t19
15\t6\t22\t14
15\t8\t16\t11.5
16\t0\t16.5\t11.5
16\t3\t27.5\t16.5
16\t4\t21.5\t15.5
16\t6\t17.5\t12
16\t8\t16\t10.5
17\t0\t12\t9.5
17\t3\t25\t19
17\t4\t21\t14.5
17\t6\t16.5\t12.5
17\t8\t14.5\t11.5
18\t0\t15.5\t11
18\t3\t30.5\t20.5
18\t4\t24.5\t17.5
18\t6\t19.5\t15
18\t8\t16.5\t11
19\t0\t12.5\t10.5
19\t1\t29\t21
19\t2\t27\t19
19\t3\t22.5\t16.5
19\t4\t22\t14
19\t6\t16\t11.5
19\t8\t13.5\t9.5
20\t0\t12\t9.5
20\t3\t18.5\t19
20\t4\t17\t17
20\t6\t16.5\t11.5
20\t8\t11\t9
21\t0\t13\t13.5
21\t3\t23\t23.5
21\t4\t21.5\t21.5
21\t6\t18.5\t17
21\t8\t15.5\t15"""

DEMO = """subjectID\tage_cat\tsex_code\tweight_cat\tdiagnosis\tALT_U_L\tAST_U_L\tHCT\tALB_g_L\tSCr_umol_L\tapache_ii\tSOFA\turine_24h_mL
1\t1\t1\t1\tAcute severe pancreatitis\t37\t52\t0.27\t34.5\t31.6\t18\t20\t66
2\t0\t1\t3\tAcute severe pancreatitis\t22\t83\t0.50\t27\t72.7\t18\t12\t153
3\t1\t1\t0\tAcute severe pancreatitis\t4\t9\t0.37\t32.7\t73\t19\t12\t89
4\t1\t1\t0\tAcute severe pancreatitis\t27\t32\t0.47\t41.7\t60.8\t22\t15\t25
5\t0\t1\t3\tAcute severe pancreatitis\t55\t93\t0.29\t30\t151.5\t5\t16\t870
6\t1\t0\t1\tAcute severe pancreatitis\t158\t193\t0.49\t45.8\t68.9\t10\t15\t50
7\t0\t1\t3\tSepsis\t37\t71\t0.26\t26.9\t163.2\t9\t17\t0
8\t1\t1\t1\tAcute severe pancreatitis\t125\t76\t0.25\t24.2\t82.5\t14\t15\t1720
9\t1\t0\t2\tSepsis\t32\t57\t0.22\t32.8\t47.6\t18\t16\t15
10\t0\t1\t1\tSepsis\t22\t66\t0.38\t27.7\t139.7\t17\t15\t510
11\t1\t1\t2\tAcute severe pancreatitis\t106\t168\t0.36\t26.7\t78.4\t21\t16\t15
12\t0\t1\t0\tAcute severe pancreatitis\t8\t19\t0.29\t28.4\t27.7\t15\t15\t505
13\t1\t1\t3\tAcute pancreatitis\t23\t43\t0.21\t26\t287\t28\t14\t0
14\t0\t1\t3\tAcute severe pancreatitis\t24\t34\t0.25\t28.5\t190.3\t15\t15\t20
15\t1\t0\t1\tAcute severe pancreatitis\t3\t15\t0.3\t30.4\t75\t5\t13\t110
16\t1\t0\t2\tAcute severe pancreatitis\t17\t39\t0.26\t29.9\t128.4\t18\t13\t75
17\t1\t0\t1\tAcute severe pancreatitis\t21\t58\t0.26\t25.5\t94.3\t19\t17\t10
18\t1\t0\t1\tAcute severe pancreatitis\t11\t21\t0.27\t32.1\t47.7\t20\t13\t10
19\t1\t1\t2\tAcute severe pancreatitis\t53\t41\t0.19\t23.8\t192.4\t12\t8\t5
20\t0\t1\t2\tAcute severe pancreatitis\t26\t29\t0.3\t27.9\t203.7\t12\t8\t170
21\t0\t1\t2\tAcute severe pancreatitis\t19\t17\t0.28\t30.4\t53.2\t27\t18\t120"""

CRRT = """subjectID\tcrrt_modality\tultrafiltration_cat\teffluent_cat\tcrrt_days_cat
1\t0\t2\t2\t1
2\t0\t1\t0\t0
3\t0\t1\t2\t1
4\t1\t1\t2\t1
5\t0\t2\t0\t3
6\t0\t1\t2\t2
7\t1\t0\t0\t3
8\t0\t0\t1\t3
9\t1\t1\t0\t0
10\t0\t0\t1\t2
11\t0\t0\t1\t3
12\t1\t2\t2\t0
13\t0\t2\t0\t0
14\t0\t0\t0\t1
15\t0\t0\t1\t2
16\t0\t2\t1\t3
17\t0\t1\t1\t3
18\t0\t0\t1\t3
19\t0\t2\t1\t1
20\t1\t1\t1\t0
21\t1\t0\t0\t1"""

out = os.path.dirname(os.path.abspath(__file__))
for name, txt in (("Ceftazidime_concentration.csv", CAZ), ("Avibactam_concentration.csv", AVI),
                  ("Demographic_data.csv", DEMO), ("CRRT_parameters_data.csv", CRRT)):
    rows = [l.split("\t") for l in txt.strip().split("\n")]
    with open(os.path.join(out, name), "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh, lineterminator="\n").writerows(rows)
    print(f"{name}: {len(rows)-1} rows")
