SELECT '>'||feature_id || ' ' ||replace(uniquename, ':', '__')
          ||regexp_replace(regexp_replace(residues, E'\\*\$', ''), E'(.{1,60})', E'\n\\1', 'g')
  FROM feature
 WHERE type_id = (
       SELECT cvterm.cvterm_id
         FROM cvterm
         JOIN cv USING (cv_id)
        WHERE cv.name = 'sequence'
          AND cvterm.name = 'polypeptide'
        )
   AND organism_id = %s

