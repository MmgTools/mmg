;;; File: emacs-format-file

;; Put the following into your ~/.emacs:
;;(load-file "PATH_TO_FILE/emacs-format-file.el")
;;(eval-after-load "emacs-format-file"
;;  '(progn
;;     (require 'emacs-format-file)
;;     (add-hook 'emacs-lisp-mode-hook 'indent-file-when-save)
;;    ))


(provide 'emacs-format-file)

;; c/c++ style
(require 'cc-mode)

;; no brace indentation
(setq-default c-default-style "linux"
              c-basic-offset 4)

(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(tab-stop-list (quote (4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100 104 108 112 116 120))))

;; no tabulation
(setq-default indent-tabs-mode nil)

;; for cmake
(eval-after-load "cmake-mode"
  '(progn
     (require 'cmake-mode)
     (setq auto-mode-alist
           (append '(("CMakeLists\\.txt\\'" . cmake-mode)
                     ("\\.cmake\\'" . cmake-mode))
                   auto-mode-alist))
     (setq cmake-tab-width 4)
     ))

;; print tab and whitespace
(require 'whitespace)
(setq whitespace-style '(face empty tabs lines-tail trailing))
(global-whitespace-mode t)

(customize-variable (quote tab-stop-list))
