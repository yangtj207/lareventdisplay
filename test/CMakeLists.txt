include(CetTest)
cet_enable_asserts()

# Tests to ensure well-formed dictionaries for the event-display
# classes (should be able to create a null pointer in ROOT).

function(evd_dictionary_test CLASS)
  cet_test(${CLASS}_t HANDBUILT
    TEST_EXEC root.exe
    TEST_ARGS -l -q -b -e "evd::${CLASS} const* x = nullptr;")
endfunction()

evd_dictionary_test(TWQProjectionView)
evd_dictionary_test(TWQMultiTPCProjectionView)
evd_dictionary_test(Display3DView)
evd_dictionary_test(Ortho3DView)
evd_dictionary_test(Ortho3DPad)
evd_dictionary_test(CalorView)
evd_dictionary_test(DrawingPad)
