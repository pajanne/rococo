SELECT cvterm.name as cvterm_name, featureprop.value as featureprop_value 
  FROM featureprop, cvterm 
 WHERE type_id = cvterm.cvterm_id 
   AND (cvterm.name = 'colour' OR cvterm.name = 'gene' OR cvterm.name = 'EC_number')
   AND featureprop.feature_id = %s
