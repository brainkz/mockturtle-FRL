Window DFS
[(a_4_s1),(b_7_s3),(c_9_s3),(d_10_m2)]
    # expanded path 
    (d_10_m2 d_9_m2)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_10_m2 d_9_m2)]
    # expanded path 
    (d_10_m2 d_9_m2 d_8_m2)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_10_m2 d_9_m2 d_8_m2)]
    # Found a complete 4-path (nphases + 1)
    (d_10_m2 d_9_m2 d_8_m2 d_7_m2)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_10_m2 d_9_m2 d_8_m2 d_7_m2)]
    # Moved into
    (d_9_m2 d_8_m2 d_7_m2 d_6_m2)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_9_m2 d_8_m2 d_7_m2 d_6_m2)]
    # Moved into two paths (due to merger m2)
    (d_8_m2 d_7_m2 d_6_m2 d_5_m2 m2_5_α)
    (d_8_m2 d_7_m2 d_6_m2 d_5_m2 m2_5_s2 s2_5_m1)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_8_m2 d_7_m2 d_6_m2 d_5_m2 m2_5_α),(d_8_m2 d_7_m2 d_6_m2 d_5_m2 m2_5_s2 s2_5_m1)]
    # Moved into
    (d_7_m2 d_6_m2 d_5_m2 m2_5_s2 s2_5_m1 s2_4_m1)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_8_m2 d_7_m2 d_6_m2 d_5_m2 m2_5_α),(d_7_m2 d_6_m2 d_5_m2 m2_5_s2 s2_5_m1 s2_4_m1)]
    # Moved into two paths (due to merger m1)
    (d_6_m2 d_5_m2 m2_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β)
    # Sigma is already clocked - path is invalid (d_6_m2 d_5_m2 m2_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_σ)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_8_m2 d_7_m2 d_6_m2 d_5_m2 m2_5_α),(d_6_m2 d_5_m2 m2_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β)]
    (d_5_m2 m2_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β s1_2_β)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_8_m2 d_7_m2 d_6_m2 d_5_m2 m2_5_α),(d_5_m2 m2_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β s1_2_β)]
    # Bumped into beta, do nothing
[(a_4_s1),(b_7_s3),(c_9_s3),(d_8_m2 d_7_m2 d_6_m2 d_5_m2 m2_5_α)]
    # Moved into
    (d_7_m2 d_6_m2 d_5_m2 m2_5_α m2_4_α)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_7_m2 d_6_m2 d_5_m2 m2_5_α m2_4_α)]
    # Moved into
    (d_6_m2 d_5_m2 m2_5_α m2_4_α m2_3_α)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_6_m2 d_5_m2 m2_5_α m2_4_α m2_3_α)]
    # Moved into
    (d_5_m2 m2_5_α m2_4_α m2_3_α m2_2_α)
[(a_4_s1),(b_7_s3),(c_9_s3),(d_5_m2 m2_5_α m2_4_α m2_3_α m2_2_α)]
    # Moved into
    (m2_5_α m2_4_α m2_3_α m2_2_α m2_1_α)
[(a_4_s1),(b_7_s3),(c_9_s3),(m2_5_α m2_4_α m2_3_α m2_2_α m2_1_α)]
    # Bumped into alpha, do nothing
[(a_4_s1),(b_7_s3),(c_9_s3)]
[(a_4_s1),(b_7_s3),(c_9_s3 c_8_s3)]
[(a_4_s1),(b_7_s3),(c_9_s3 c_8_s3 c_7_s3)]
    # Found complete 4-path 
    (c_9_s3 c_8_s3 c_7_s3 c_6_s3 s3_6_s2)
[(a_4_s1),(b_7_s3),(c_9_s3 c_8_s3 c_7_s3 c_6_s3 s3_6_s2)]
    # Moved into
    (c_8_s3 c_7_s3 c_6_s3 s3_6_s2 s3_5_s2 s2_5_m1)
[(a_4_s1),(b_7_s3),(c_8_s3 c_7_s3 c_6_s3 s3_6_s2 s3_5_s2 s2_5_m1)]
    # Moved into
    (c_7_s3 c_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1)
[(a_4_s1),(b_7_s3),(c_7_s3 c_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1)]
    # Moved into
    (c_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β)
    # PATH INVALID σ==m1 (c_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_σ)
[(a_4_s1),(b_7_s3),(c_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β)]
    # Moved into
    (s3_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β s1_2_β)
[(a_4_s1),(b_7_s3),(s3_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β s1_2_β)]
    # Bumped into beta, do nothing
[(a_4_s1),(b_7_s3)]
    # extended path
    (b_7_s3 b_6_s3 s3_6_s2)
[(a_4_s1),(b_7_s3 b_6_s3 s3_6_s2)]
    # extended path
    (b_7_s3 b_6_s3 s3_6_s2 s3_5_s2 s2_5_m1)
[(a_4_s1),(b_7_s3 b_6_s3 s3_6_s2 s3_5_s2 s2_5_m1)]
    # Found complete 4-path
    (b_7_s3 b_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1)
[(a_4_s1),(b_7_s3 b_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1)]
    # Moved into
    (b_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β)
    # PATH INVALID (b_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_σ)
[(a_4_s1),(b_6_s3 s3_6_s2 s3_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β)]
    # Moved into
    # ALREADY SEEN (s3_5_s2 s2_5_m1 s2_4_m1 s2_3_m1 m1_3_s1 s1_3_β s1_2_β)
[(a_4_s1)]
    # extended path
    (a_4_s1 a_3_s1 s1_3_β)
[(a_4_s1 a_3_s1 s1_3_β)]
    # extended path
    (a_4_s1 a_3_s1 s1_3_β s1_2_β)
[(a_4_s1 a_3_s1 s1_3_β s1_2_β)]
    # bumped into beta, do nothing
[]