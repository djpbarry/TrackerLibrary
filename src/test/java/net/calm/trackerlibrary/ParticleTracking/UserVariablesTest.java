package net.calm.trackerlibrary.ParticleTracking;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;

class UserVariablesTest {

    @AfterEach
    void resetGlobalState() {
        UserVariables.setInstance(null);
    }

    @Test
    void getInstanceIsASingletonAndCanBeReplaced() {
        UserVariables a = UserVariables.getInstance();
        UserVariables b = UserVariables.getInstance();
        assertSame(a, b);

        UserVariables fresh = new UserVariables();
        UserVariables.setInstance(fresh);
        assertSame(fresh, UserVariables.getInstance());
        assertNotSame(a, UserVariables.getInstance());
    }

    @Test
    void instanceIsolationKeepsSeparateTrajectoriesFromBleeding() {
        UserVariables clean = new UserVariables();
        clean.setTimeRes(7.5);
        UserVariables.setInstance(clean);

        assertEquals(7.5, UserVariables.getInstance().getTimeRes(), 1e-12);
        // A second instance does not inherit the mutation.
        assertEquals(1.0, new UserVariables().getTimeRes(), 1e-12);
    }
}
