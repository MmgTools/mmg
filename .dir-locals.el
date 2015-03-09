;; Set emacs indentation

((nil . (
         (fill-column              . 80 )
         (mode                     . whitespace)
         (indent-tabs-mode         . nil)
         (tab-width                . 4 )
         (show-trailing-whitespace . t  )
         )
      )

 (c-mode . (
            (indent-tabs-mode         . nil)
            (c-file-style             . "linux")
            (c-basic-offset           . 4  )))

 (cmake-mode . (
                (cmake-file-style            . "linux")
                (cmake-tab-width          . 4  )))
 )
