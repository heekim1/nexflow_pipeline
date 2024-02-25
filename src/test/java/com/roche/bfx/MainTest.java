package com.roche.bfx;

import com.roche.bfx.error.BfxError;
import com.roche.bfx.sbx.consensus.SomeComponent;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertNotNull;

class MainTest {

    /**
     * Stub to test the visibility of the common component.
     */
    @Test
    void viewCommon() {
        BfxError bfxError = new BfxError() {
        };
        assertNotNull(bfxError);
    }

    /**
     * Stub to test visibility of other components.
     */
    @Test
    void viewComponent() {
        SomeComponent someComponent = new SomeComponent();
        assertNotNull(someComponent);
    }

}
